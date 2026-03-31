"""
@author: cunyue
@file: test_git.py
@time: 2026/3/31
@description: Tests for git information collection
"""

from subprocess import CalledProcessError, TimeoutExpired
from unittest.mock import MagicMock, patch

import pytest

from swanlab.sdk.internal.run.system.environment import git
from swanlab.sdk.typings.run.system import GitSnapshot


class TestParseGitUrl:
    """Tests for git.parse_git_url()"""

    @pytest.mark.parametrize(
        "url,expected",
        [
            ("git@github.com:SwanHubX/SwanLab.git", "https://github.com/SwanHubX/SwanLab.git"),
            ("git@gitlab.com:user/project.git", "https://gitlab.com/user/project.git"),
            ("git@gitlab.com:group/subgroup/project.git", "https://gitlab.com/group/subgroup/project.git"),
            ("git@github.com:user/repo", "https://github.com/user/repo"),
            ("https://github.com/user/repo.git", "https://github.com/user/repo.git"),
            ("http://github.com/user/repo.git", "http://github.com/user/repo.git"),
            ("git://github.com/user/repo.git", "git://github.com/user/repo.git"),
            ("", ""),
        ],
    )
    def test_parse_git_url(self, url, expected):
        """测试各种 URL 格式转换"""
        assert git.parse_git_url(url) == expected


class TestGetRemoteUrl:
    """Tests for git.get_remote_url()"""

    @pytest.mark.parametrize(
        "stdout,expected",
        [
            ("https://github.com/user/repo.git\n", "https://github.com/user/repo.git"),
            ("git@github.com:user/repo.git\n", "https://github.com/user/repo.git"),
        ],
    )
    def test_get_remote_url_success(self, stdout, expected):
        """成功获取远程仓库地址"""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout=stdout)
            result = git.get_remote_url()

        assert result == expected

    @pytest.mark.parametrize(
        "exception",
        [
            CalledProcessError(1, "git"),
            TimeoutExpired("git", 5),
            FileNotFoundError("git not found"),
        ],
    )
    def test_get_remote_url_returns_none_on_error(self, exception):
        """各种错误情况下返回 None"""
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = exception
            result = git.get_remote_url()

        assert result is None


class TestGetBranch:
    """Tests for git.get_branch()"""

    @pytest.mark.parametrize(
        "stdout,expected",
        [
            ("main\n", "main"),
            ("feature/new-feature\n", "feature/new-feature"),
            ("", ""),
        ],
    )
    def test_get_branch_success(self, stdout, expected):
        """成功获取分支名"""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout=stdout)
            result = git.get_branch()

        assert result == expected

    @pytest.mark.parametrize(
        "exception",
        [
            CalledProcessError(1, "git"),
            TimeoutExpired("git", 5),
            FileNotFoundError("git not found"),
        ],
    )
    def test_get_branch_returns_none_on_error(self, exception):
        """各种错误情况下返回 None"""
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = exception
            result = git.get_branch()

        assert result is None


class TestGetCommit:
    """Tests for git.get_commit()"""

    def test_get_commit_success(self):
        """成功获取当前提交 hash"""
        with (
            patch("swanlab.sdk.internal.run.system.environment.git.get_branch") as mock_branch,
            patch("subprocess.run") as mock_run,
        ):
            mock_branch.return_value = "main"
            mock_run.return_value = MagicMock(returncode=0, stdout="abc123def456\n")
            result = git.get_commit()

        assert result == "abc123def456"

    @pytest.mark.parametrize("branch_value", [None, ""])
    def test_get_commit_no_branch(self, branch_value):
        """没有分支时返回 None"""
        with patch("swanlab.sdk.internal.run.system.environment.git.get_branch") as mock_branch:
            mock_branch.return_value = branch_value
            result = git.get_commit()

        assert result is None

    @pytest.mark.parametrize(
        "exception",
        [
            CalledProcessError(1, "git"),
            TimeoutExpired("git", 5),
        ],
    )
    def test_get_commit_returns_none_on_error(self, exception):
        """各种错误情况下返回 None"""
        with (
            patch("swanlab.sdk.internal.run.system.environment.git.get_branch") as mock_branch,
            patch("subprocess.run") as mock_run,
        ):
            mock_branch.return_value = "main"
            mock_run.side_effect = exception
            result = git.get_commit()

        assert result is None


class TestGitGet:
    """Tests for git.get()"""

    @pytest.mark.parametrize(
        "url,branch,commit",
        [
            ("https://github.com/user/repo.git", "main", "abc123"),
            (None, None, None),
            ("https://github.com/user/repo.git", None, None),
        ],
    )
    def test_get_returns_git_snapshot(self, url, branch, commit):
        """get() 返回 GitSnapshot 实例"""
        with (
            patch("swanlab.sdk.internal.run.system.environment.git.get_remote_url", return_value=url),
            patch("swanlab.sdk.internal.run.system.environment.git.get_branch", return_value=branch),
            patch("swanlab.sdk.internal.run.system.environment.git.get_commit", return_value=commit),
        ):
            result = git.get()

        assert isinstance(result, GitSnapshot)
        assert result.remote_url == url
        assert result.branch == branch
        assert result.commit == commit
