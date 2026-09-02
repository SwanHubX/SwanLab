"""Temporary compatibility APIs.

These APIs are intentionally isolated from the regular SDK surface and may be
changed or removed without compatibility or performance guarantees.
"""

import warnings
from typing import Any, Dict, Mapping

from swanlab.api import Api
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.run import Run, fmt, get_run

from .guard import with_cmd_lock, with_run


def _get_missing_data(run: Run, data: Mapping[str, Any], step: int) -> Dict[str, Any]:
    """Return only scalar metrics that have no point at ``step`` on the server."""
    flattened_data = fmt.flatten_dict(data)
    if not flattened_data:
        return {}

    query_result = Api().run(run.path.strip("/")).metrics(
        keys=list(flattened_data),
        range_query={"start": step, "end": step},
    )
    points_by_key = {item["key"]: item["metrics"] for item in query_result["list"]}
    missing_keys = set(flattened_data) - set(points_by_key)
    if missing_keys:
        raise RuntimeError(f"SwanLab did not return metrics for {sorted(missing_keys)!r}; refusing to backfill")

    missing_data: Dict[str, Any] = {}
    for key, value in flattened_data.items():
        if any(point["step"] == step for point in points_by_key[key]):
            console.warning(f"Skip backfill: metric {key!r} at step={step} already exists on SwanLab.")
        else:
            missing_data[key] = value
    return missing_data


def _allow_backfill(run: Run, data: Mapping[str, Any], step: int) -> None:
    """Allow this one historical step through the Python core's local guard."""
    metrics = getattr(run._ctx.core, "_metrics", None)
    if metrics is None:
        return

    for key in fmt.flatten_dict(data):
        metric = metrics.get(key)
        if metric is not None:
            # Only relax the guard, never make it stricter. Removing this exact
            # step from the local bitmap also permits a retry in the same run.
            metric.min_step = min(metric.min_step, step - 1)
            metric.steps.discard(step)


@with_cmd_lock
@with_run("log_backfill")
def log_backfill(data: Mapping[str, Any], step: int) -> None:
    """Temporarily log historical data that Resume would otherwise skip.

    .. warning::

       This is a temporary compatibility API. Its behavior, availability, and
       performance are not guaranteed and it may be changed or removed.

    Only scalar metrics are supported. The active run must be an online run
    created with ``resume=True``. The server is queried before uploading:
    existing points are skipped with a warning, while missing points are
    uploaded. ``step`` must be an explicit, non-negative integer.
    """
    warnings.warn(
        "swanlab.log_backfill() is a temporary compatibility API; its behavior, "
        "availability, and performance are not guaranteed.",
        FutureWarning,
        stacklevel=2,
    )

    run = get_run()
    if isinstance(data, Mapping) and isinstance(step, int) and step >= 0:
        if run._ctx.config.settings.run.resume == "never":
            raise RuntimeError("swanlab.log_backfill() requires a run initialized with resume=True.")
        missing_data = _get_missing_data(run, data, step)
        if not missing_data:
            return
        _allow_backfill(run, missing_data, step)
        run.log(missing_data, step=step)
        return
    run.log(data, step=step)
