from datetime import datetime, timezone

import pytest

from swanlab.utils.time import TIMESTAMP_S_MAX, TIMESTAMP_S_MIN, parse_cli_timestamp_ns, parse_resp_timestamp_s


class TestParseRespTimestampS:
    def test_valid_variants(self):
        """测试秒/毫秒时间戳与 ISO 8601 字符串的归一化"""
        assert parse_resp_timestamp_s(1722470400) == 1722470400
        assert parse_resp_timestamp_s("1722470400") == 1722470400
        assert parse_resp_timestamp_s(1722470400000) == 1722470400
        assert parse_resp_timestamp_s("1722470400000") == 1722470400
        assert parse_resp_timestamp_s("2024-08-01T00:00:00Z") == 1722470400
        assert parse_resp_timestamp_s("2024-08-01T00:00:00.123Z") == 1722470400
        assert parse_resp_timestamp_s("2024-08-01T00:00:00+00:00") == 1722470400
        assert parse_resp_timestamp_s("2024-08-01T08:00:00+08:00") == 1722470400
        # 无时区时按 UTC 解析
        assert parse_resp_timestamp_s("2024-08-01T00:00:00") == 1722470400

    def test_valid_boundaries(self):
        """测试合法区间边界值"""
        assert parse_resp_timestamp_s(TIMESTAMP_S_MIN) == TIMESTAMP_S_MIN
        assert parse_resp_timestamp_s(TIMESTAMP_S_MAX) == TIMESTAMP_S_MAX

    @pytest.mark.parametrize(
        "value",
        [None, "", "   ", "not-a-date", 0, -100, "0", "-100", 100, "999999999", "1970-01-01T00:00:00Z"],
    )
    def test_invalid_raises(self, value):
        """缺失 / 空值 / 无法解析 / 非正数 / 越界必须抛错，禁止回退默认值导致慢查询"""
        with pytest.raises(ValueError):
            parse_resp_timestamp_s(value)


class TestParseCliTimestampNs:
    def test_accepts_iso_and_unix(self):
        """ISO 8601 字符串与秒/毫秒时间戳都归一化为同一个纳秒时间戳"""
        expected = 1722470400 * 1_000_000_000

        assert parse_cli_timestamp_ns("2024-08-01T00:00:00Z") == expected
        assert parse_cli_timestamp_ns(1722470400) == expected
        assert parse_cli_timestamp_ns(1722470400000) == expected
        assert parse_cli_timestamp_ns("1722470400") == expected

    def test_keeps_subsecond_precision(self):
        """亚秒输入必须精确到纳秒，否则闭区间端点会被判到范围外。"""
        base = 1722470400 * 1_000_000_000

        assert parse_cli_timestamp_ns("2024-08-01T00:00:00.000001Z") == base + 1_000
        assert parse_cli_timestamp_ns("2024-08-01T00:00:00.000002Z") == base + 2_000
        assert parse_cli_timestamp_ns(datetime(2024, 8, 1, 0, 0, 0, 1, tzinfo=timezone.utc)) == base + 1_000
        assert parse_cli_timestamp_ns(1722470400123) == base + 123_000_000
        assert parse_cli_timestamp_ns("1722470400123") == base + 123_000_000
        assert parse_cli_timestamp_ns(1722470400.001) == base + 1_000_000

    def test_treats_naive_datetime_as_local(self):
        """不带时区的输入按系统本地时区解释，与终端输入墙上时间的习惯一致"""
        naive = datetime(2024, 8, 1, 8, 30, 0)

        assert parse_cli_timestamp_ns(naive) == int(round(naive.timestamp() * 1_000_000_000))
        assert parse_cli_timestamp_ns("2024-08-01 08:30:00") == int(round(naive.timestamp() * 1_000_000_000))

    @pytest.mark.parametrize("value", ["", "  ", "not-a-time", -1, 0, float("nan")])
    def test_rejects_invalid_values(self, value):
        with pytest.raises(ValueError):
            parse_cli_timestamp_ns(value)

    @pytest.mark.parametrize(
        "value",
        [20261020, "20261020", "1970-01-01T00:00:01Z", datetime(1960, 1, 1, tzinfo=timezone.utc)],
    )
    def test_rejects_out_of_range(self, value):
        """漏写分隔符的日期（20261020）不能被当作秒静默放行，否则过滤条件等于失效。"""
        with pytest.raises(ValueError, match="out of range"):
            parse_cli_timestamp_ns(value)

    def test_accepts_lower_bound(self):
        assert parse_cli_timestamp_ns(TIMESTAMP_S_MIN) == TIMESTAMP_S_MIN * 1_000_000_000
        assert parse_cli_timestamp_ns(TIMESTAMP_S_MAX) == TIMESTAMP_S_MAX * 1_000_000_000

    def test_rejects_bool(self):
        with pytest.raises(TypeError):
            parse_cli_timestamp_ns(True)
