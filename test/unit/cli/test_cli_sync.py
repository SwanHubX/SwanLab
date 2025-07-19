"""
@author: cunyue
@file: test_cli_sync.py
@time: 2025/7/19 20:01
@description: 测试 sync 功能
"""

import nanoid
import pytest
from click.testing import CliRunner

from swanlab.cli.commands.sync import sync


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def temp_dirs(tmp_path):
    """创建两个临时目录"""
    dir1 = tmp_path / "dir1"
    dir2 = tmp_path / "dir2"
    dir1.mkdir()
    dir2.mkdir()
    return str(dir1), str(dir2)


@pytest.fixture
def temp_id():
    """创建一个临时 ID"""
    return nanoid.generate(alphabet="0123456789abcdefghijklmnopqrstuvwxyz", size=21)


# noinspection PyTypeChecker
class TestClicSyncResume:
    """
    测试 resume 相关规则
    - 如果使用了 `--resume` 选项，则只能同步一个目录。
    - 使用了 `--resume` 选项，则不能同时使用 `--id`、`--workspace` 或 `--project` 选项。
    - 如果使用了 `--id` 选项，则只能同步一个目录。
    """

    @staticmethod
    def test_resume_multiple_paths(runner, temp_dirs):
        """resume + 多个路径 - 报错"""
        result = runner.invoke(sync, [temp_dirs[0], temp_dirs[1], "--resume"])
        assert result.exit_code == 2
        assert "only be used with a single directory path" in result.output

    @staticmethod
    def test_id_multiple_paths(runner, temp_dirs, temp_id):
        """id + 多个路径 - 报错"""
        result = runner.invoke(sync, [temp_dirs[0], temp_dirs[1], "--id", temp_id])
        assert result.exit_code == 2
        assert "only be used with a single directory path" in result.output

    @staticmethod
    def test_resume_and_id_multiple_paths(runner, temp_dirs, temp_id):
        """resume + id - 报错"""
        result = runner.invoke(sync, [temp_dirs[0], temp_dirs[1], "--resume", "--id", temp_id])
        assert result.exit_code == 2
        assert "The --resume option cannot be used with --id, --workspace, or --project options." in result.output
