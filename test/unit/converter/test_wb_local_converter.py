"""
@file: test_wb_local_converter.py
@description: Unit tests for WandB local converter
"""

import json
import os
import tempfile
from unittest.mock import Mock, patch, MagicMock
import pytest

from swanlab.converter.wb.wb_local_converter import WandbLocalConverter


class TestWandbLocalConverterInit:
    """Test WandbLocalConverter initialization"""

    def test_init_with_defaults(self):
        converter = WandbLocalConverter()
        assert converter.project is None
        assert converter.workspace is None
        assert converter.mode == "cloud"
        assert converter.tags is None
        assert converter.logdir is None

    def test_init_with_custom_params(self):
        converter = WandbLocalConverter(
            project="test-project",
            workspace="test-workspace",
            mode="local",
            tags=["tag1", "tag2"],
            logdir="/tmp/logs"
        )
        assert converter.project == "test-project"
        assert converter.workspace == "test-workspace"
        assert converter.mode == "local"
        assert converter.tags == ["tag1", "tag2"]
        assert converter.logdir == "/tmp/logs"


class TestFindRunDirs:
    """Test _find_run_dirs method"""

    def test_find_run_dirs_with_specific_run(self, tmp_path):
        root_dir = tmp_path / "wandb"
        root_dir.mkdir()
        run_dir = root_dir / "run-12345"
        run_dir.mkdir()

        converter = WandbLocalConverter()
        result = converter._find_run_dirs(str(root_dir), "run-12345")

        assert len(result) == 1
        assert str(run_dir) in result[0]

    def test_find_run_dirs_with_pattern(self, tmp_path):
        root_dir = tmp_path / "wandb"
        root_dir.mkdir()
        (root_dir / "run-12345").mkdir()
        (root_dir / "run-67890").mkdir()
        (root_dir / "offline-run-11111").mkdir()

        converter = WandbLocalConverter()
        result = converter._find_run_dirs(str(root_dir))

        assert len(result) == 3

    def test_find_run_dirs_no_matches(self, tmp_path):
        root_dir = tmp_path / "wandb"
        root_dir.mkdir()

        converter = WandbLocalConverter()
        result = converter._find_run_dirs(str(root_dir))

        assert len(result) == 0


class TestUnpackKeyValueJsonList:
    """Test _unpack_key_value_json_list method"""

    def test_unpack_valid_key_value_list(self):
        converter = WandbLocalConverter()
        items = [
            {"key": "learning_rate", "value_json": "0.001"},
            {"key": "batch_size", "value_json": "32"}
        ]

        result = converter._unpack_key_value_json_list(items)

        assert result == {"learning_rate": 0.001, "batch_size": 32}

    def test_unpack_nested_key(self):
        converter = WandbLocalConverter()
        items = [
            {"nested_key": ["model", "layers"], "value_json": "12"}
        ]

        result = converter._unpack_key_value_json_list(items)

        assert result == {"model/layers": 12}

    def test_unpack_invalid_json(self):
        converter = WandbLocalConverter()
        items = [
            {"key": "valid", "value_json": "123"},
            {"key": "invalid", "value_json": "{invalid json}"}
        ]

        result = converter._unpack_key_value_json_list(items)

        assert "valid" in result
        assert result["valid"] == 123
        assert "invalid" not in result

    def test_unpack_non_list_returns_as_is(self):
        converter = WandbLocalConverter()
        items = {"not": "a list"}

        result = converter._unpack_key_value_json_list(items)

        assert result == items

    def test_unpack_empty_list(self):
        converter = WandbLocalConverter()
        result = converter._unpack_key_value_json_list([])
        assert result == []


class TestRunMethod:
    """Test run method"""

    @patch('swanlab.converter.wb.wb_local_converter.WandbLocalConverter._parse_run')
    @patch('swanlab.converter.wb.wb_local_converter.WandbLocalConverter._find_run_dirs')
    def test_run_with_no_runs_found(self, mock_find, mock_parse):
        mock_find.return_value = []

        converter = WandbLocalConverter()
        converter.run("/fake/path")

        mock_find.assert_called_once()
        mock_parse.assert_not_called()

    @patch('swanlab.converter.wb.wb_local_converter.WandbLocalConverter._parse_run')
    @patch('swanlab.converter.wb.wb_local_converter.WandbLocalConverter._find_run_dirs')
    def test_run_with_multiple_runs(self, mock_find, mock_parse):
        mock_find.return_value = ["/path/run1", "/path/run2"]

        converter = WandbLocalConverter()
        converter.run("/fake/path")

        assert mock_parse.call_count == 2

    @patch('swanlab.converter.wb.wb_local_converter.WandbLocalConverter._parse_run')
    @patch('swanlab.converter.wb.wb_local_converter.WandbLocalConverter._find_run_dirs')
    def test_run_handles_parse_exception(self, mock_find, mock_parse):
        mock_find.return_value = ["/path/run1", "/path/run2"]
        mock_parse.side_effect = [Exception("Parse error"), None]

        converter = WandbLocalConverter()
        converter.run("/fake/path")

        assert mock_parse.call_count == 2
