from unittest.mock import MagicMock, patch

import pytest

from swanlab.sdk.internal.probe_python.hardware_vendor.accelerator.huawei import (
    AscendNPU,
)

ASCEND_950_MAPPING = """\
NPU ID  Slot ID  Chip ID  Chip Phy-ID  Chip Name
0       1        0        0            Ascend950PR
1       2        0        1            Ascend950PR
"""

LEGACY_MAPPING = """\
NPU ID  Chip ID  Chip Logic ID  Chip Name
0       0        0              Ascend910B
0       1        MCU            Ascend910B
"""


@pytest.mark.parametrize(
    ("output", "expected"),
    [
        (
            ASCEND_950_MAPPING,
            {
                "0": {"0": {"id": "0", "name": "Ascend950PR", "chip_query": False}},
                "1": {"0": {"id": "1", "name": "Ascend950PR", "chip_query": False}},
            },
        ),
        (
            LEGACY_MAPPING,
            {"0": {"0": {"id": "0", "name": "Ascend910B", "chip_query": True}}},
        ),
    ],
)
def test_map_npu_raw_supports_mapping_formats(output, expected):
    with patch("subprocess.run", return_value=MagicMock(stdout=output)):
        assert AscendNPU._map_npu_raw() == expected


def test_ascend_950_queries_without_chip_id():
    outputs = {
        "usages": """\
HBM Capacity(MB)               : 131072
HBM Usage Rate(%)              : 4
Aicore Usage Rate(%)           : 37
""",
        "temp": "Temperature (C)                : 62\n",
        "power": "NPU Real-time Power(W)         : 203.3\n",
    }

    def run(command, **_kwargs):
        return MagicMock(stdout=outputs[command[3]])

    collector = object.__new__(AscendNPU)
    with patch("subprocess.run", side_effect=run) as mock_run:
        usage = collector._collect_usage("0", "0", "0-0", chip_query=False)
        temp = collector._collect_temp("0", "0", "0-0", chip_query=False)
        power = collector._collect_power("0", "0", "0-0", chip_query=False)

    assert usage == [
        ("npu.0-0.pct", 37.0),
        ("npu.0-0.mem.pct", 4.0),
        ("npu.0-0.mem.value", 5242.88),
    ]
    assert temp == ("npu.0-0.temp", 62.0)
    assert power == ("npu.0-0.power", 203.3)
    assert all("-c" not in call.args[0] for call in mock_run.call_args_list)


def test_legacy_query_includes_chip_id():
    with patch("subprocess.run", return_value=MagicMock(stdout="")) as mock_run:
        AscendNPU._query_info("usages", "0", "1")

    assert mock_run.call_args.args[0] == ["npu-smi", "info", "-t", "usages", "-i", "0", "-c", "1"]
