import pytest
from click import BadParameter

from swanlab.cli.api.helper import validate_filter_query


def test_inline_json():
    result = validate_filter_query('[{"key":"state","type":"STABLE","op":"EQ","value":["RUNNING"]}]')
    assert len(result) == 1
    assert result[0]["key"] == "state"


def test_from_file(tmp_path):
    p = tmp_path / "f.json"
    p.write_text('[{"key":"name","type":"STABLE","op":"CONTAIN","value":["test"]}]')
    result = validate_filter_query(str(p))
    assert len(result) == 1
    assert result[0]["op"] == "CONTAIN"


def test_empty_raises():
    with pytest.raises(BadParameter, match="must not be empty"):
        validate_filter_query("  ")


def test_invalid_json_raises():
    with pytest.raises(BadParameter, match="neither a valid file path nor valid JSON"):
        validate_filter_query("not json")


def test_non_array_raises():
    with pytest.raises(BadParameter, match="JSON array"):
        validate_filter_query('{"key": "state"}')
