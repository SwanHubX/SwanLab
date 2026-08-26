#!/usr/bin/env bash
# SwanLab 多平台 release 构建脚本。
#
# 产物：1 sdist + 1 any 兜底 wheel + 6 平台 wheel（内嵌交叉编译的 Go core）。
# 交叉编译与平台 tag 写入由 hatch_build.py 构建钩子完成，
# 本脚本负责切环境变量逐平台构建，以及产物完整性校验。
#
# 用法：
#   bash scripts/build_release.sh [VERSION]      # 更新版本号并构建全部产物，随后校验
#   bash scripts/build_release.sh                # 使用 swanlab/package.json 现有版本
#   bash scripts/build_release.sh --verify-only  # 仅校验已有 dist/，不构建
set -euo pipefail
cd "$(dirname "$0")/.."

PLATFORMS=(
  "manylinux_2_17_x86_64.manylinux2014_x86_64"   # linux/amd64
  "manylinux_2_17_aarch64.manylinux2014_aarch64" # linux/arm64
  "win_amd64"                                    # windows/amd64
  "win_arm64"                                    # windows/arm64
  "macosx_12_0_x86_64"                           # darwin/amd64
  "macosx_12_0_arm64"                            # darwin/arm64
)

VERIFY_ONLY=false
VERSION=""
for arg in "$@"; do
  case "$arg" in
    --verify-only) VERIFY_ONLY=true ;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) VERSION="$arg" ;;
  esac
done

if [ "$VERIFY_ONLY" = false ]; then
  # 0. 版本号（沿用 make build 的写回逻辑）
  if [ -n "$VERSION" ]; then
    python3 -c "import json; data=json.load(open('swanlab/package.json')); data['version']='$VERSION'; json.dump(data,open('swanlab/package.json','w'),indent=2)"
    echo "==> Updated swanlab/package.json version to $VERSION"
  fi

  # 1. 清理旧产物
  rm -rf dist swanlab/bin

  # 2. sdist + any 兜底 wheel（不设 SWANLAB_BUILD_PLATFORM，钩子零动作）
  echo "==> Building sdist + fallback any wheel"
  uv build --sdist
  uv build --wheel

  # 3. 逐平台交叉编译（钩子内完成：go build + 原生平台 tag + artifacts 注入）
  for TAG in "${PLATFORMS[@]}"; do
    echo "==> Building platform wheel: $TAG"
    SWANLAB_BUILD_PLATFORM="$TAG" uv build --wheel
    rm -rf swanlab/bin # 清理，避免污染下一平台
  done
fi

# 4. 校验
echo "==> twine check"
uvx twine check dist/*

echo "==> 结构校验"
python3 - <<'PY'
import pathlib
import re
import subprocess
import sys
import tarfile
import tempfile
import zipfile

DIST = pathlib.Path("dist")

# tag -> (二进制文件名, file 输出需命中的关键词，任一命中即可)
PLATFORMS = {
    "manylinux_2_17_x86_64.manylinux2014_x86_64": ("swanlab-core", ["elf", "x86-64"]),
    "manylinux_2_17_aarch64.manylinux2014_aarch64": ("swanlab-core", ["elf", "aarch64"]),
    "win_amd64": ("swanlab-core.exe", ["pe32+", "x86-64", "x86_64"]),
    "win_arm64": ("swanlab-core.exe", ["pe32+", "aarch64", "arm64"]),
    "macosx_12_0_x86_64": ("swanlab-core", ["mach-o", "x86_64"]),
    "macosx_12_0_arm64": ("swanlab-core", ["mach-o", "arm64"]),
}

failures = []


def check(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg)
    if not cond:
        failures.append(msg)


def wheel_version(name):
    m = re.fullmatch(r"swanlab-(.+)-py3-none-(.+)\.whl", name)
    return (m.group(1), m.group(2)) if m else None


wheels = sorted(DIST.glob("swanlab-*-py3-none-*.whl"))
sdists = sorted(DIST.glob("swanlab-*.tar.gz"))
parsed = [wheel_version(w.name) for w in wheels]
check(all(p is not None for p in parsed), "wheel 文件名均可解析")
check(len(sdists) == 1, f"恰好 1 个 sdist（实际 {len(sdists)}）")

# 版本一致性：所有产物（含 sdist）使用同一版本串
versions = {p[0] for p in parsed if p} | {re.fullmatch(r"swanlab-(.+)\.tar\.gz", s.name).group(1) for s in sdists}
check(len(versions) == 1, f"全部产物版本一致（实际 {versions}）")

# 产物齐全：any 兜底 + 6 平台 wheel，无多余文件
tags = {p[1] for p in parsed if p}
check(tags == set(PLATFORMS) | {"any"}, f"wheel tag 齐全（多出: {sorted(tags - set(PLATFORMS) - {'any'})}, 缺少: {sorted(set(PLATFORMS) - tags)}）")

for whl in wheels:
    tag = wheel_version(whl.name)[1]
    with zipfile.ZipFile(whl) as zf:
        names = zf.namelist()
        has_binary = any(n.startswith("swanlab/bin/") for n in names)

        if tag == "any":
            check(not has_binary, "any wheel 不含二进制")
            continue

        exe = PLATFORMS[tag][0]
        entry = f"swanlab/bin/{exe}"
        if entry not in names:
            check(False, f"{tag}: 缺少 {entry}")
            continue

        mode = zf.getinfo(entry).external_attr >> 16
        check(mode & 0o111 != 0, f"{tag}: 二进制含可执行位")

        with tempfile.TemporaryDirectory() as td:
            extracted = zf.extract(entry, td)
            desc = subprocess.run(["file", "-b", str(extracted)], capture_output=True, text=True).stdout.lower()
        check(any(k in desc for k in PLATFORMS[tag][1]), f"{tag}: file 架构匹配（{desc.strip()}）")

        if tag.startswith("manylinux"):
            # glibc tag 双写：WHEEL 元数据需逐条展开为多行 Tag:
            meta_name = next(n for n in names if n.endswith(".dist-info/WHEEL"))
            meta = zf.read(meta_name).decode()
            lines = [ln.split(": ", 1)[1] for ln in meta.splitlines() if ln.startswith("Tag:")]
            expected = {f"py3-none-{t}" for t in tag.split(".")}
            check(expected <= set(lines), f"{tag}: WHEEL 元数据多行 Tag 展开（实际 {lines}）")

# sdist：不含 core/ 源码与二进制，但必须携带构建钩子（保证从 sdist 重建 wheel 可行）
with tarfile.open(sdists[0]) as tf:
    members = tf.getnames()
check(any(m.count("/") == 1 and m.endswith("/hatch_build.py") for m in members), "sdist 携带 hatch_build.py")
check(not any(m.count("/") == 1 and m.split("/")[1] == "core" for m in members), "sdist 不含 core/ 源码")
check(not any("/swanlab/bin/" in m for m in members), "sdist 不含二进制")

if failures:
    print(f"\n{len(failures)} 项校验失败", file=sys.stderr)
    sys.exit(1)
print("\nAll checks passed.")
PY

echo "==> Done"
