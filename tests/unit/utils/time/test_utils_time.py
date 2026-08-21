import pytest

from swanlab.utils.time import TIMESTAMP_S_MAX, TIMESTAMP_S_MIN, parse_timestamp_s


class TestParseTimestampS:
    def test_valid_variants(self):
        """测试秒/毫秒时间戳与 ISO 8601 字符串的归一化"""
        assert parse_timestamp_s(1722470400) == 1722470400
        assert parse_timestamp_s("1722470400") == 1722470400
        assert parse_timestamp_s(1722470400000) == 1722470400
        assert parse_timestamp_s("1722470400000") == 1722470400
        assert parse_timestamp_s("2024-08-01T00:00:00Z") == 1722470400
        assert parse_timestamp_s("2024-08-01T00:00:00.123Z") == 1722470400
        assert parse_timestamp_s("2024-08-01T00:00:00+00:00") == 1722470400
        assert parse_timestamp_s("2024-08-01T08:00:00+08:00") == 1722470400
        # 无时区时按 UTC 解析
        assert parse_timestamp_s("2024-08-01T00:00:00") == 1722470400

    def test_valid_boundaries(self):
        """测试合法区间边界值"""
        assert parse_timestamp_s(TIMESTAMP_S_MIN) == TIMESTAMP_S_MIN
        assert parse_timestamp_s(TIMESTAMP_S_MAX) == TIMESTAMP_S_MAX

    @pytest.mark.parametrize(
        "value",
        [None, "", "   ", "not-a-date", 0, -100, "0", "-100", 100, "999999999", "1970-01-01T00:00:00Z"],
    )
    def test_invalid_raises(self, value):
        """缺失 / 空值 / 无法解析 / 非正数 / 越界必须抛错，禁止回退默认值导致慢查询"""
        with pytest.raises(ValueError):
            parse_timestamp_s(value)
