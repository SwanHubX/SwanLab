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

    # Base command for protoc
    command = [
        sys.executable,
        "-m",
        "grpc_tools.protoc",
        f"-I{proto_dir}",
        f"--python_out={output_dir}",
        f"--grpc_python_out={output_dir}",
    ]

    command.extend([str(p.relative_to(proto_dir)) for p in proto_files])

    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error generating Python protos:\n{result.stderr}")
        sys.exit(1)
    else:
        print("Python protos generated successfully.")

    for root, dirs, files in os.walk(output_dir):
        root_path = Path(root)
        if not (root_path / "__init__.py").exists():
            (root_path / "__init__.py").touch()


def generate_go_protos(proto_dir: Path, output_dir: Path):
    """
    Placeholder and guide for Go proto generation within the 'core' module.
    """
    print("\n[Go Generation Guide (Monorepo: swanlab-core)]")
    print("To generate Go protos for the 'core' module, run:")

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
