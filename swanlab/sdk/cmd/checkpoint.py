"""
@author: Zirui Cai
@file: checkpoint.py
@time: 2026/11/25
@description: swanlab.checkpoint 顶层函数：登记一个 checkpoint 到平台
"""

from typing import Dict, List, Optional, cast

from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.core_python.api.checkpoint import register_checkpoint
from swanlab.sdk.internal.core_python.core import CorePython
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.run import get_run, has_run

__all__ = ["checkpoint"]

# sentinel：区分「未传 step」和「显式传 None（留空）」
_UNSET = object()


def checkpoint(
    uri: str,
    *,
    run_id: Optional[str] = None,
    step=_UNSET,  # type: ignore[assignment]  # .pyi 中声明为 Optional[int]
    name: Optional[str] = None,
    description: Optional[str] = None,
    metadata: Optional[Dict] = None,
    tags: Optional[List[str]] = None,
    username: Optional[str] = None,
    project: Optional[str] = None,
) -> str:
    """Register a checkpoint with the SwanLab platform.

    Records the URI of a checkpoint (local path, S3, Hugging Face Hub URL, etc.)
    along with metadata about the training state at that point.

    :param uri: Path or URI pointing to the checkpoint. Any scheme is accepted
        (local path, ``s3://``, ``hf://``, etc.). The platform stores this value
        as-is without accessing or validating the target.
    :param run_id: Source experiment run ID (cuid). Defaults to the current active
        run when one exists. Required if no run is active.
    :param step: Training step at which the checkpoint was saved. Defaults to the
        current global step of the active run. Pass ``None`` explicitly to leave
        the step unset in the platform.
    :param name: Human-readable name for the checkpoint. If omitted, the platform
        derives one from the last segment of ``uri``.
    :param description: Optional free-text note.
    :param metadata: Optional dict of arbitrary key/value pairs (e.g. optimizer
        settings, learning rate) stored alongside the checkpoint record.
    :param tags: Optional list of tag names, e.g. ``["best", "ready-for-eval"]``.
    :param username: Username of the project owner. Auto-resolved from the active
        run when one exists. Required if no run is active.
    :param project: Project name. Auto-resolved from the active run when one exists.
        Required if no run is active.
    :return: The platform-assigned checkpoint ``cuid``, or an empty string when
        registration is skipped (e.g. in non-online mode without an active run).

    Examples:

        Register inside a training loop (step, username, and project resolved automatically):

        >>> import swanlab
        >>> run = swanlab.init(project="my_project")
        >>> for step in range(1000):
        ...     run.log({"loss": loss})
        ...     if step % 100 == 0:
        ...         torch.save(model.state_dict(), f"ckpt-{step}.pt")
        ...         swanlab.checkpoint(f"./ckpt-{step}.pt")

        Register with explicit run_id (e.g. in a separate script):

        >>> cuid = swanlab.checkpoint(
        ...     uri="s3://my-bucket/run-a/checkpoint-500",
        ...     run_id="clxi2v9abc0001...",
        ...     username="my-org",
        ...     project="my-project",
        ...     step=500,
        ...     name="val-best",
        ...     tags=["best"],
        ... )
    """
    active_run = get_run() if has_run() else None

    # ── 非 online 模式：无法登记（在访问 core._ctx 之前检查）──────────────────
    if active_run is not None and active_run.mode != "online":
        console.warning(
            f"swanlab.checkpoint() is only available in online mode (current mode: {active_run.mode}). "
            "This checkpoint will not be registered."
        )
        return ""

    # ── 无活跃 run 且无 client（无法发请求）────────────────────────────────────
    if active_run is None and not client.exists():
        console.warning(
            "swanlab.checkpoint() requires an active online run or an initialized client. "
            "This checkpoint will not be registered."
        )
        return ""

    # ── 解析 core context（online 且有活跃 run 时）────────────────────────────
    core_ctx = None
    if active_run is not None:
        core = active_run._ctx.core
        assert isinstance(core, CorePython), "checkpoint() requires an online CorePython instance"
        core_ctx = core._ctx

    # ── 解析 experiment_id ────────────────────────────────────────────────────
    if run_id is not None:
        experiment_id = run_id
    elif core_ctx is not None:
        experiment_id = core_ctx.experiment_id
    else:
        console.warning(
            "swanlab.checkpoint() requires a run_id. "
            "Either call swanlab.init() first or pass run_id= explicitly. "
            "This checkpoint will not be registered."
        )
        return ""

    # ── 解析 username / project ───────────────────────────────────────────────
    resolved_username = username or (core_ctx.username if core_ctx is not None else None)
    resolved_project = project or (core_ctx.project if core_ctx is not None else None)
    if not resolved_username or not resolved_project:
        console.warning(
            "swanlab.checkpoint() requires username and project. "
            "Either call swanlab.init() first or pass username= and project= explicitly. "
            "This checkpoint will not be registered."
        )
        return ""

    # ── 解析 step ─────────────────────────────────────────────────────────────
    if step is _UNSET:
        gs = active_run._ctx._global_step if active_run is not None else -1
        resolved_step: Optional[int] = gs if gs >= 0 else None
    else:
        resolved_step = cast(Optional[int], step)

    return register_checkpoint(
        username=resolved_username,
        project=resolved_project,
        experiment_id=experiment_id,
        uri=uri,
        step=resolved_step,
        name=name,
        description=description,
        metadata=metadata,
        tags=tags,
    )
