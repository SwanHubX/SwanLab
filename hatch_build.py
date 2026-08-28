"""SwanLab wheel 自定义构建钩子（hatchling 插件）。

由 ``SWANLAB_BUILD_PLATFORM`` 环境变量控制构建行为：

- 未设置：不执行任何操作，产出纯 ``py3-none-any`` wheel。
- 设置为平台 tag（如 ``manylinux_2_17_aarch64.manylinux2014_aarch64``、
  ``win_arm64``、``macosx_12_0_arm64``）：交叉编译 Go core 内嵌进 wheel，
  并将该平台 tag 写入 wheel 元数据。
"""

import os
import pathlib
import re
import shutil
import sys

from hatchling.builders.hooks.plugin.interface import BuildHookInterface

_BUILD_PLATFORM_ENV = "SWANLAB_BUILD_PLATFORM"

# 支持的平台 tag 形态示例（GOOS 由前缀判定，GOARCH 由结尾后缀判定）
_SUPPORTED_TAG_EXAMPLES = (
    "manylinux_2_17_x86_64.manylinux2014_x86_64",
    "manylinux_2_17_aarch64.manylinux2014_aarch64",
    "win_amd64",
    "win_arm64",
    "macosx_12_0_x86_64",
    "macosx_12_0_arm64",
)


class CustomBuildHook(BuildHookInterface):
    """按目标平台编译 Go core，并产出对应平台 tag 的 wheel。"""

    def initialize(self, version, build_data):
        if self.target_name == "wheel":
            self._prepare_wheel(build_data)

    def _prepare_wheel(self, build_data):
        platform_tag = os.getenv(_BUILD_PLATFORM_ENV)
        if not platform_tag:
            return

        goos, goarch = self._parse_target_platform(platform_tag)
        output = self._build_core(goos, goarch)

        build_data["tag"] = f"py3-none-{platform_tag}"
        # hatchling 统一使用正斜杠路径，Windows 构建环境下亦然
        build_data["artifacts"].append(output.as_posix())

    def _parse_target_platform(self, platform_tag):
        """从平台 tag 解析 (GOOS, GOARCH)，无法解析时中止构建。"""
        tag = platform_tag.strip().lower()
        if tag.startswith(("manylinux_", "musllinux_", "linux_")):
            goos = "linux"
        elif tag.startswith(("win_", "win-")):
            goos = "windows"
        elif tag.startswith(("macosx_", "macosx-")):
            goos = "darwin"
        else:
            self.app.abort(
                f"Unrecognized {_BUILD_PLATFORM_ENV}={platform_tag!r}: cannot determine OS prefix. "
                f"Supported tags: {', '.join(_SUPPORTED_TAG_EXAMPLES)}",
            )
            raise AssertionError("unreachable")

        arch_match = re.search(r"(?:-|_)(aarch64|arm64|x86_64|amd64)$", tag)
        if not arch_match:
            self.app.abort(
                f"Cannot parse target architecture from {_BUILD_PLATFORM_ENV}={platform_tag!r} "
                "(expected suffix _amd64 / _x86_64 / _arm64 / _aarch64)",
            )
            raise AssertionError("unreachable")

        goarch = {"amd64": "amd64", "x86_64": "amd64", "arm64": "arm64", "aarch64": "arm64"}[arch_match.group(1)]
        return goos, goarch

    def _build_core(self, goos, goarch):
        """编译 Go core 到 swanlab/bin/，返回产物路径（相对仓库根）。"""
        go = shutil.which("go")
        if not go:
            self.app.abort(
                f"{_BUILD_PLATFORM_ENV} is set but no Go toolchain found. "
                "Install Go to build platform wheels: https://go.dev/doc/install",
            )
            raise AssertionError("unreachable")

        output = pathlib.Path("swanlab", "bin", "swanlab-core")
        if goos == "windows":
            output = output.with_suffix(".exe")

        self.app.display_waiting(f"Building swanlab-core ({goos}/{goarch})...")
        try:
            # 惰性导入：sdist 不携带 core/ 源码，仅平台构建时才需要该模块
            sys.path.insert(0, str(pathlib.Path(__file__).parent))
            from core import hatch as hatch_core
        except ImportError:
            self.app.abort(
                f"{_BUILD_PLATFORM_ENV} is set but core/hatch.py is missing. "
                "Platform wheels must be built from a git checkout; sdist does not include Go sources.",
            )
            raise AssertionError("unreachable")

        hatch_core.build_core(
            go_binary=pathlib.Path(go),
            output_path=output,
            target_system=goos,
            target_arch=goarch,
        )
        return output
