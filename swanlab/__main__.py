"""
@description: Package entry point for ``python -m swanlab``.
"""

from swanlab.cli import cli


def main() -> None:
    cli(prog_name="python -m swanlab")


if __name__ == "__main__":
    main()
