"""
@author: Nexisato
@file: test_errors.py
@time: 2026/10/27
@description: CoreGoError / ReasonCode / reason_for_grpc 的单元测试。
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage, reason_for_grpc

REPO_ROOT = Path(__file__).resolve().parents[5]


def test_error_str_includes_stage_reason_and_grpc_name():
    error = CoreGoError(
        Stage.RPC,
        ReasonCode.UNIMPLEMENTED,
        "UpsertScalars is not implemented yet",
        grpc_name="UNIMPLEMENTED",
        grpc_code=12,
        detail="server details",
    )

    assert str(error) == "[RPC/UNIMPLEMENTED] UpsertScalars is not implemented yet (grpc=UNIMPLEMENTED)"
    assert isinstance(error, Exception)
    assert error.args == ("UpsertScalars is not implemented yet",)


def test_error_str_without_grpc_name():
    error = CoreGoError(Stage.SPAWN, ReasonCode.BINARY_NOT_FOUND, "binary missing")

    assert str(error) == "[SPAWN/BINARY_NOT_FOUND] binary missing"


def test_error_to_dict_is_json_friendly():
    error = CoreGoError(
        Stage.PORT_FILE,
        ReasonCode.OWNERSHIP_MISMATCH,
        "pid mismatch",
        detail="stderr tail",
    )

    assert error.to_dict() == {
        "stage": "PORT_FILE",
        "reason": "OWNERSHIP_MISMATCH",
        "message": "pid mismatch",
        "grpc_name": None,
        "grpc_code": None,
        "detail": "stderr tail",
    }


def test_error_keeps_original_cause():
    cause = RuntimeError("boom")
    error = CoreGoError(Stage.CHANNEL, ReasonCode.GRPC_UNAVAILABLE, "connect failed", cause=cause)

    assert error.__cause__ is cause


@pytest.mark.parametrize(
    "grpc_name, expected",
    [
        ("UNIMPLEMENTED", ReasonCode.UNIMPLEMENTED),
        ("UNAVAILABLE", ReasonCode.GRPC_UNAVAILABLE),
        ("DEADLINE_EXCEEDED", ReasonCode.DEADLINE_EXCEEDED),
        ("FAILED_PRECONDITION", ReasonCode.FAILED_PRECONDITION),
        ("NOT_FOUND", ReasonCode.SESSION_NOT_FOUND),
        ("RESOURCE_EXHAUSTED", ReasonCode.RESOURCE_EXHAUSTED),
        ("INVALID_ARGUMENT", ReasonCode.INVALID_ARGUMENT),
        ("INTERNAL", ReasonCode.GRPC_ERROR),
        ("", ReasonCode.GRPC_ERROR),
    ],
)
def test_reason_for_grpc_mapping(grpc_name, expected):
    assert reason_for_grpc(grpc_name) is expected


def test_reason_for_grpc_is_case_insensitive():
    assert reason_for_grpc(" unimplemented ") is ReasonCode.UNIMPLEMENTED


def test_import_core_go_does_not_load_grpc():
    """core_go 顶层导入不触发 grpc；用子进程避免测试自身提前加载。"""
    code = "import sys; import swanlab.sdk.internal.core_go; print('grpc' in sys.modules)"
    env = dict(os.environ, PYTHONPATH=str(REPO_ROOT))

    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )

    assert result.stdout.strip() == "False", result.stderr
