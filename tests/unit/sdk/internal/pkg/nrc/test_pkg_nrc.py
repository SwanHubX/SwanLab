import netrc as _netrc
import sys
from pathlib import Path

import pytest

from swanlab.sdk.internal.pkg import nrc


class TestFmt:
    @pytest.mark.parametrize(
        "i, expected",
        [
            ("https://example.com/path?a=1", "https://example.com"),
            ("https://example.com/a/b/c?key=val#hash", "https://example.com"),
            ("http://example.com:8080/api", "http://example.com:8080"),
            ("example.com:3000/test", "https://example.com:3000"),
            ("example.com", "https://example.com"),
            ("http://example.com", "http://example.com"),
            ("https://example.com", "https://example.com"),
            ("http://example.com:9090", "http://example.com:9090"),
            ("example.com:5000", "https://example.com:5000"),
            ("  https://example.com/path  ", "https://example.com"),
            ("://", "://"),
        ],
    )
    def test_fmt(self, i, expected):
        assert nrc.fmt(i) == expected

    @pytest.mark.parametrize("i", ["", "   "])
    def test_empty_raises(self, i):
        with pytest.raises(ValueError):
            nrc.fmt(i)


class TestRead:
    def test_file_not_exists(self, tmp_path):
        assert nrc.read(tmp_path / ".netrc") is None

    def test_empty_file(self, tmp_path):
        p = tmp_path / ".netrc"
        p.write_text("")
        assert nrc.read(p) is None

    def test_normal(self, tmp_path):
        p = tmp_path / ".netrc"
        p.write_text("machine api.swanlab.cn login swanlab.cn password test_key")
        assert nrc.read(p) == ("test_key", "https://api.swanlab.cn", "https://swanlab.cn")

    def test_with_scheme(self, tmp_path):
        p = tmp_path / ".netrc"
        p.write_text("machine https://api.swanlab.cn login https://swanlab.cn password test_key")
        assert nrc.read(p) == ("test_key", "https://api.swanlab.cn", "https://swanlab.cn")

    def test_with_port(self, tmp_path):
        p = tmp_path / ".netrc"
        p.write_text("machine private.com:8080 login private.com password my_key")
        assert nrc.read(p) == ("my_key", "https://private.com:8080", "https://private.com")

    def test_legacy_api_suffix(self, tmp_path):
        """向下兼容：旧版本写入了带 /api 后缀的 machine，应自动修复"""
        p = tmp_path / ".netrc"
        p.write_text("machine api.swanlab.cn/api login swanlab.cn password old_key")
        assert nrc.read(p) == ("old_key", "https://api.swanlab.cn", "https://swanlab.cn")
        parsed = _netrc.netrc(p)
        assert "api.swanlab.cn/api" not in parsed.hosts
        assert "api.swanlab.cn" in parsed.hosts

    def test_corrupted_file(self, tmp_path):
        p = tmp_path / ".netrc"
        p.write_text("invalid content without machine keyword")
        with pytest.raises(IOError):
            nrc.read(p)

    def test_directory_conflict(self, tmp_path):
        p = tmp_path / ".netrc"
        p.mkdir()
        with pytest.raises(IOError):
            nrc.read(p)


class TestWrite:
    def test_normal(self, tmp_path):
        p = tmp_path / ".netrc"
        nrc.write(p, "api.swanlab.cn", "swanlab.cn", "test_key")
        assert p.exists()
        auth = _netrc.netrc(p).authenticators("api.swanlab.cn")
        assert auth is not None
        assert (auth[0], auth[2]) == ("swanlab.cn", "test_key")

    def test_permissions(self, tmp_path):
        p = tmp_path / ".netrc"
        nrc.write(p, "api.swanlab.cn", "swanlab.cn", "test_key")
        if sys.platform != "win32":
            assert (p.stat().st_mode & 0o777) == 0o600

    def test_sso_overwrites_previous(self, tmp_path):
        """全局单点登录：后写入的必须清空之前不同环境的凭证"""
        p = tmp_path / ".netrc"
        nrc.write(p, "api.swanlab.cn", "swanlab.cn", "key1")
        nrc.write(p, "private.com", "private.com", "key2")
        parsed = _netrc.netrc(p)
        assert "api.swanlab.cn" not in parsed.hosts
        assert len(parsed.hosts) == 1
        auth = parsed.authenticators("private.com")
        assert auth is not None
        assert (auth[0], auth[2]) == ("private.com", "key2")

    def test_skip_if_same(self, tmp_path):
        """凭证相同时跳过写入"""
        p = tmp_path / ".netrc"
        nrc.write(p, "api.swanlab.cn", "swanlab.cn", "test_key")
        mtime_before = p.stat().st_mtime
        nrc.write(p, "api.swanlab.cn", "swanlab.cn", "test_key")
        assert p.stat().st_mtime == mtime_before

    def test_recovery_on_corrupted(self, tmp_path):
        """文件损坏时能重置并成功写入"""
        p = tmp_path / ".netrc"
        p.write_text("broken file contents!!")
        nrc.write(p, "api.swanlab.cn", "swanlab.cn", "new_key")
        auth = _netrc.netrc(p).authenticators("api.swanlab.cn")
        assert auth is not None
        assert (auth[0], auth[2]) == ("swanlab.cn", "new_key")

    def test_directory_conflict(self, tmp_path):
        p = tmp_path / ".netrc"
        p.mkdir()
        with pytest.raises(IOError, match="is a directory"):
            nrc.write(p, "api.swanlab.cn", "swanlab.cn", "key")

    def test_invalid_characters(self):
        with pytest.raises(ValueError, match="Invalid characters"):
            nrc.write(Path("/tmp/.netrc"), "api.swanlab.cn", "swanlab.cn", "bad key")

    def test_write_then_read_roundtrip(self, tmp_path):
        """写入后读取应保持数据一致"""
        p = tmp_path / ".netrc"
        nrc.write(p, "api.swanlab.cn", "swanlab.cn", "round_trip_key")
        assert nrc.read(p) == ("round_trip_key", "https://api.swanlab.cn", "https://swanlab.cn")
