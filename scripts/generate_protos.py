"""
@author: cunyue
@file: generate_protos.py
@time: 2026/3/11 18:36
@description: 生成 SwanLab ProtoBuf 文件

对于Go部分，你需要安装：

1. brew install protobuf
2. go install google.golang.org/protobuf/cmd/protoc-gen-go@v1.36.11
3. go install google.golang.org/grpc/cmd/protoc-gen-go-grpc@v1.6.1

"""

import os
import re
import subprocess
import sys
from pathlib import Path


def generate_python_protos(proto_dir: Path, output_dir: Path):
    """
    Generate Python code to swanlab/proto.
    """
    print(f"Generating Python protos from {proto_dir} to {output_dir}...")

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all .proto files
    proto_files = list(proto_dir.glob("**/*.proto"))
    if not proto_files:
        print("No .proto files found.")
        return

    # 1. 执行生成
    command = [
        sys.executable,
        "-m",
        "grpc_tools.protoc",
        f"-I{proto_dir}",
        f"--python_out={output_dir}",
        f"--grpc_python_out={output_dir}",
        f"--pyi_out={output_dir}",
    ]

    command.extend([str(p.relative_to(proto_dir)) for p in proto_files])

    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error generating Python protos:\n{result.stderr}")
        sys.exit(1)
    else:
        print("Python protos generated successfully.")

    # 2. 修复导入路径 (SwanLab 特供逻辑)
    # 我们需要将 'from swanlab.xxx.v1 import ...'
    # 替换为 'from swanlab.proto.swanlab.xxx.v1 import ...'
    print("Fixing import paths for internal SDK usage...")

    # 匹配模式：针对你展示的生成的代码格式
    # 寻找以 from swanlab 开头的导入行
    import_re = re.compile(r"^from swanlab\.(.*) import (.*) as (.*)$", re.MULTILINE)

    for py_file in output_dir.rglob("*.py"):
        if py_file.name == "__init__.py":
            continue

        with open(py_file, "r", encoding="utf-8") as f:
            content = f.read()

        # 核心修复：插入 .proto 路径层级
        # 结果会变成: from swanlab.proto.swanlab.run.v1 import run_pb2 as ...
        fixed_content = import_re.sub(r"from swanlab.proto.swanlab.\1 import \2 as \3", content)

        # 兼容性修复：处理可能存在的 'import xxx_pb2' 相对引用
        fixed_content = re.sub(r"^import (.*_pb2)", r"from . import \1", fixed_content, flags=re.MULTILINE)

        with open(py_file, "w", encoding="utf-8") as f:
            f.write(fixed_content)

    # 3. 补全 __init__.py (确保整个 proto 树都是可导入的 package)
    for root, _, _ in os.walk(output_dir):
        init_file = Path(root) / "__init__.py"
        if not init_file.exists():
            init_file.touch()


def generate_go_protos(proto_dir: Path, output_dir: Path):
    """
    Placeholder and guide for Go proto generation within the 'core' module.
    """
    # In a monorepo, we typically want the generated files to be within the module directory.
    # --go_opt=paths=source_relative ensures the directory structure mirrors the source.
    go_cmd = [
        "protoc",
        f"-I{proto_dir}",
        f"--go_out={output_dir}",
        "--go_opt=paths=source_relative",
        f"--go-grpc_out={output_dir}",
        "--go-grpc_opt=paths=source_relative",
    ]

    proto_files = [str(p.relative_to(proto_dir)) for p in proto_dir.glob("**/*.proto")]
    go_cmd.extend(proto_files)

    try:
        result = subprocess.run(go_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error generating Go protos:\n{result.stderr}")
            # 提示：如果报错 command not found，说明没装 protoc-gen-go
            if "not found" in result.stderr or "executable file not found" in result.stderr:
                print(
                    "\nTip: Install Go plugins via:\n  go install google.golang.org/protobuf/cmd/protoc-gen-go@latest"
                )
        else:
            print("Go protos generated successfully.")
    except FileNotFoundError:
        print("Error: 'protoc' command not found. Please install protoc first.")


if __name__ == "__main__":
    project_root = Path(__file__).parent.parent

    # Source: protos/
    proto_source = project_root / "protos"

    # Python: swanlab/proto/
    python_output = project_root / "swanlab" / "proto"

    # Go: core/proto/ (Part of the swanlab-core module)
    go_output = project_root / "core" / "proto"

    # Ensure core/proto exists for the guide/placeholder
    go_output.mkdir(parents=True, exist_ok=True)

    # Run generation
    generate_python_protos(proto_source, python_output)
    generate_go_protos(proto_source, go_output)
