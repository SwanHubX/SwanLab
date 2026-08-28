"""swanlab-core Go 模块编译封装。

供两处调用：

- ``hatch_build.py`` 构建钩子（发布/平台 wheel 构建）
- 直接执行 ``python3 core/hatch.py``（``make core-bin``，本机构建）

职责单一：以 ``core/`` 为工作目录执行 ``go build``，将产物输出到仓库根的
``swanlab/bin/swanlab-core(.exe)``，并注入版本号与 commit。

注：``CGO_ENABLED=0`` 静态编译下，Go 内部链接器写入的 ELF ``.gnu.version`` /
``.gnu.version_r`` section 会导致 auditwheel 崩溃。仅当未来引入 manylinux
容器认证（cibuildwheel 路线）时，才需要在构建流程中用 objcopy 移除这两个
section；当前自声明 manylinux tag 的构建不涉及 auditwheel，无需处理。
"""

from __future__ import annotations

import json
import os
import pathlib
import subprocess

_REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
_CORE_DIR = _REPO_ROOT / "core"
_ENTRY_PACKAGE = "./cmd/swanlab-core"


def build_core(
    go_binary: pathlib.Path,
    output_path: pathlib.PurePath,
    target_system: str | None = None,
    target_arch: str | None = None,
) -> None:
    """编译 swanlab-core。

    Args:
        go_binary: go 可执行文件路径，必须存在。
        output_path: 产物路径，相对仓库根（如 ``swanlab/bin/swanlab-core``）。
        target_system: 目标 GOOS，``None`` 表示使用本机平台。
        target_arch: 目标 GOARCH，``None`` 表示使用本机平台。
    """
    version = _package_version()
    commit = _git_commit()
    # -s -w 移除符号表与 DWARF 调试信息；-X 将版本信息注入 main 包变量
    ldflags = f"-s -w -X main.version={version} -X main.commit={commit}"
    # go build 以 core/ 为工作目录，因此输出路径需要相对 core/ 前移一级
    output = pathlib.Path("..") / output_path

    subprocess.check_call(
        [
            str(go_binary),
            "build",
            "-trimpath",
            f"-ldflags={ldflags}",
            "-o",
            str(output),
            _ENTRY_PACKAGE,
        ],
        cwd=str(_CORE_DIR),
        env=_go_env(target_system, target_arch),
    )

    if not _is_windows_target(target_system):
        # 显式设置可执行位：wheel 以 zip 外部属性记录权限，pip 解压时还原
        os.chmod(_REPO_ROOT / output_path, 0o755)


def _go_env(target_system: str | None, target_arch: str | None) -> dict[str, str]:
    env = os.environ.copy()
    if target_system:
        env["GOOS"] = target_system
    if target_arch:
        env["GOARCH"] = target_arch
    # 纯 Go 无 cgo：静态编译，保证全平台交叉编译无宿主依赖
    env["CGO_ENABLED"] = "0"
    return env


def _is_windows_target(target_system: str | None) -> bool:
    if target_system:
        return target_system == "windows"
    return os.name == "nt"


def _package_version() -> str:
    """从 swanlab/package.json 读取版本号，与 SDK 版本保持同源。"""
    package_json = _REPO_ROOT / "swanlab" / "package.json"
    version = json.loads(package_json.read_text(encoding="utf-8"))["version"]
    return str(version)


def _git_commit() -> str:
    """返回当前 commit SHA；不在 git 仓库中或 git 不可用时返回 ``unknown``。"""
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(_REPO_ROOT), text=True).strip()
    except Exception:
        return "unknown"


if __name__ == "__main__":
    # 本机构建入口（make core-bin），复用与发布构建一致的编译参数
    exe = "swanlab-core.exe" if os.name == "nt" else "swanlab-core"
    go = pathlib.Path(os.environ.get("GO", "go"))
    build_core(go_binary=go, output_path=pathlib.PurePath("swanlab", "bin", exe))
    print(f"swanlab/bin/{exe} built for host platform")
