import time

import psutil

from swanlab.data.run.metadata.hardware.disk import (
    DiskCollector,
)


class TestDiskCollector:
    collector = DiskCollector()

    def test_disk_read_speed(self):
        """
        测试硬盘信息采集
        """
        start = time.time()
        current_disk_io = psutil.disk_io_counters()
        time.sleep(1)
        result = self.collector.get_disk_read_speed(current_disk_io.read_bytes, time.time() - start)
        assert result['name'] == "MB read from disk"
        assert result['value'] >= 0
        assert result['config'].chart_name == 'Disk I/O Utilization (MB)'

    def test_disk_write_speed(self):
        """
        测试硬盘信息采集
        """
        start = time.time()
        current_disk_io = psutil.disk_io_counters()
        time.sleep(1)
        result = self.collector.get_disk_write_speed(current_disk_io.write_bytes, time.time() - start)
        assert result['name'] == "MB written to disk"
        assert result['value'] >= 0
        assert result['config'].chart_name == 'Disk I/O Utilization (MB)'

    def test_disk_usage(self):
        """
        测试硬盘信息采集
        """
        result = self.collector.get_disk_usage()
        assert result['name'] == "Disk Utilization (%)"
        assert result['value'] >= 0
        assert result['config'].chart_name == 'Disk Utilization (%)'
