"""Tests for BiMap adapter"""

import pytest

from swanlab.sdk.internal.pkg.adapter.bimap import BiMap


def test_bimap_forward_mapping():
    """Test forward mapping (external -> internal)"""
    adapter = BiMap({"a": 1, "b": 2})
    assert adapter["a"] == 1
    assert adapter["b"] == 2


def test_bimap_reverse_mapping():
    """Test reverse mapping (internal -> external)"""
    adapter = BiMap({"a": 1, "b": 2})
    assert adapter[1] == "a"
    assert adapter[2] == "b"


def test_bimap_key_error():
    """Test KeyError when mapping not found"""
    adapter = BiMap({"a": 1})
    with pytest.raises(KeyError):
        _ = adapter["nonexistent"]


def test_bimap_get_with_default():
    """Test get() method with default value"""
    adapter = BiMap({"a": 1})
    assert adapter.get("a") == 1
    assert adapter.get("nonexistent", "default") == "default"
    assert adapter.get("nonexistent") is None
