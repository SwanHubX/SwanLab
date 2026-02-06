import time

from swanlab.data.porter.mounter import Mounter


class TestMounterMetrics:
    def test_get_metrics_basic(self):
        """
        Test basic functionality of get_metrics with scalar and media summaries.
        """
        columns = [
            {"key": "loss", "type": "float", "class": "SCALAR"},
            {"key": "image", "type": "image", "class": "MEDIA"},
        ]
        summaries = {
            "scalar": [{"key": "loss", "step": 10}],
            "media": [{"key": "image", "step": 5}],
        }

        metrics = Mounter.get_metrics(columns, summaries)

        assert metrics["loss"] == ("float", "SCALAR", None, 10)
        assert metrics["image"] == ("image", "MEDIA", None, 5)

    def test_get_metrics_priority(self):
        """
        Test that scalar summaries take precedence over media summaries for the same key.
        """
        columns = [
            {"key": "mixed", "type": "float", "class": "SCALAR"},
        ]
        summaries = {
            "scalar": [{"key": "mixed", "step": 20}],
            "media": [{"key": "mixed", "step": 15}],
        }

        metrics = Mounter.get_metrics(columns, summaries)

        assert metrics["mixed"] == ("float", "SCALAR", None, 20)

    def test_get_metrics_missing_summary(self):
        """
        Test that columns with no corresponding summary get step -1.
        """
        columns = [
            {"key": "empty", "type": "float", "class": "SCALAR"},
        ]
        summaries = {
            "scalar": [],
            "media": [],
        }

        metrics = Mounter.get_metrics(columns, summaries)

        assert metrics["empty"] == ("float", "SCALAR", None, -1)

    def test_get_metrics_system_filtering(self):
        """
        Test that non-SDK system metrics are filtered out.
        """
        columns = [
            {"key": "system/cpu", "type": "float", "class": "SYSTEM"},  # Valid system key
            {"key": "custom_system", "type": "float", "class": "SYSTEM"},  # Invalid system key
            {"key": "normal", "type": "float", "class": "SCALAR"},
        ]
        summaries = {
            "scalar": [
                {"key": "system/cpu", "step": 100},
                {"key": "custom_system", "step": 50},
                {"key": "normal", "step": 1},
            ],
        }

        # Mounter uses swanlab.data.run.metadata.hardware.is_system_key which checks a predefined list.
        # Since we are running in unit tests, we rely on the actual implementation of is_system_key.
        # "system/cpu" is NOT in the default system keys list usually (it's hardware specific naming like cpu_usage),
        # but let's check what is_system_key actually checks.

        # To avoid dependency on the exact list of system keys which might change or be empty in test env,
        # we can mock is_system_key. But Mounter imports it directly.
        # So we should probably check what `is_system_key` does or mock it in the test.

        # Let's mock is_system_key for this test to be deterministic.
        from unittest.mock import patch

        with patch("swanlab.data.porter.mounter.is_system_key") as mock_is_system_key:
            # Setup mock: return True for "system/cpu", False for others
            mock_is_system_key.side_effect = lambda k: k == "system/cpu"

            metrics = Mounter.get_metrics(columns, summaries)

            assert "system/cpu" in metrics
            assert "custom_system" not in metrics
            assert "normal" in metrics

    def test_get_metrics_error_handling(self):
        """
        Test handling of columns with errors.
        """
        columns = [
            {"key": "broken", "type": "float", "class": "SCALAR", "error": {"msg": "failed"}},
        ]
        summaries = {}

        metrics = Mounter.get_metrics(columns, summaries)

        assert metrics["broken"] == ("float", "SCALAR", {"msg": "failed"}, -1)

    def test_get_metrics_large_scale(self, count=1_000_000):
        """
        Test performance of get_metrics with 1 million columns.
        """

        # Generate 1 million columns
        columns = [{"key": f"metric_{i}", "type": "float", "class": "SCALAR"} for i in range(count)]

        # Generate summaries for half of them
        # 500k scalar summaries
        summaries = {
            "scalar": [{"key": f"metric_{i}", "step": i} for i in range(0, count, 2)],
            "media": [],
        }

        start_time = time.time()

        metrics = Mounter.get_metrics(columns, summaries)

        end_time = time.time()
        duration = end_time - start_time

        # Verify correctness for a few items
        assert metrics["metric_0"][3] == 0
        assert metrics["metric_2"][3] == 2
        assert metrics["metric_1"][3] == -1  # odd numbers have no summary

        # Performance assertion: Should be very fast (e.g., < 2s on modern machine)
        # Setting a conservative limit of 10s to account for CI variance
        assert duration < 10.0, f"Performance too slow: {duration}s"
