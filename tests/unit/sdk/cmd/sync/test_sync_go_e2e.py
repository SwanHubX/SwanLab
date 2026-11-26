"""
@author: Nexisato
@file: test_sync_go_e2e.py
@time: 2026/11/18
@description: 默认 Go core 路由下 swanlab.init(offline) -> swanlab.sync() 端到端集成测试。
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import swanlab
from swanlab.sdk.cmd.login import login_raw
from swanlab.sdk.cmd.sync import sync as sync_run
from swanlab.sdk.internal import impl
from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.internal.settings.experiment import ProjectSettings, RunSettings


def _create_offline_run_dir(*, steps: int = 5, finish_run: bool = True, log_dir: Optional[Path] = None) -> Path:
    settings = Settings(log_dir=log_dir) if log_dir is not None else Settings()
    run = swanlab.init(mode="offline", project="sync-go-e2e", settings=settings)
    run_dir = Path(run._ctx.run_dir)
    for step in range(1, steps + 1):
        swanlab.log({"loss": round(1.0 / step, 4)}, step=step)
    if finish_run:
        run.finish()
    else:
        from swanlab.sdk.internal.run import clear_run

        core = run._ctx.core
        clear_run()
        core.abort()
    return run_dir


def test_sync_default_go_full_with_active_client_api_key_fallback(
    tmp_path: Path, require_core_binary, start_mock_cloud, assert_no_owned_core_process
):
    assert impl.resolve_core_sync_backend().impl == "go"

    run_dir = _create_offline_run_dir(steps=5, log_dir=tmp_path / "logs")
    with start_mock_cloud(use_request_cuid=False) as cloud:
        state = cloud.state
        url = cloud.url
        login_raw("session-only-key", host=url, print_welcome=False)
        sync_run(
            run_dir,
            settings=Settings(api_host=url, web_host=url, project=ProjectSettings(name="sync-go-e2e")),
        )

    assert_no_owned_core_process()
    assert state.exp_creates == 1
    assert state.stop_states == ["FINISHED"]
    assert sorted(state.scalar_steps_by_exp["exp-1"]) == [1, 2, 3, 4, 5]
    assert "loss" in state.series_keys


def test_sync_default_go_incremental_range_and_corrupt_tail(
    tmp_path: Path, require_core_binary, start_mock_cloud, assert_no_owned_core_process
):
    run_dir = _create_offline_run_dir(steps=20, log_dir=tmp_path / "logs")
    with start_mock_cloud(
        existing_cuids={"existing-run": "existing-run"},
        summary={"scalar": [{"key": "loss", "step": 12}]},
        use_request_cuid=False,
    ) as cloud:
        state = cloud.state
        url = cloud.url
        sync_settings = Settings(
            api_host=url,
            web_host=url,
            api_key="smoke-key",
            project=ProjectSettings(name="sync-go-e2e"),
            run=RunSettings(id="existing-run"),
        )
        # 1. 增量同步：mock summary 为 step=12，仅应补传 13..20
        sync_run(run_dir, settings=sync_settings)
        assert sorted(state.scalar_steps_by_exp["existing-run"]) == list(range(13, 21))

        # 2. 范围回放：start_step=5, end_step=10，忽略 summary 水位 12 重放 5..10
        state.scalar_steps_by_exp.clear()
        sync_run(run_dir, settings=sync_settings, start_step=5, end_step=10)
        assert sorted(state.scalar_steps_by_exp["existing-run"]) == list(range(5, 11))

        # 3. 损坏尾部容错：截断末尾 4 字节后仍可同步前缀记录
        run_files = list(run_dir.glob("run-*.swanlab"))
        assert len(run_files) == 1
        raw = run_files[0].read_bytes()
        run_files[0].write_bytes(raw[:-4])

        corrupt_dir = tmp_path / "corrupt-copy"
        corrupt_dir.mkdir(parents=True, exist_ok=True)
        for p in run_dir.iterdir():
            if p.is_file():
                (corrupt_dir / p.name).write_bytes(p.read_bytes())

        state.scalar_steps_by_exp.clear()
        state.stop_states.clear()
        sync_run(
            corrupt_dir,
            settings=Settings(
                api_host=url, web_host=url, api_key="smoke-key", project=ProjectSettings(name="sync-go-e2e")
            ),
        )
        uploaded = state.scalar_steps_by_exp["exp-1"]
        assert len(uploaded) > 0
        assert uploaded[0] == 1
        assert len(uploaded) <= 20
        # 截断掉了尾部 FinishRecord，新实验按 CRASHED 收尾
        assert state.stop_states == ["CRASHED"]

    assert_no_owned_core_process()
