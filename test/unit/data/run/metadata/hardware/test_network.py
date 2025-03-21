import time

import psutil

from swanlab.data.run.metadata.hardware.network import (
    NetworkCollector,
)


class TestNetworkCollector:
    collector = NetworkCollector()

    def test_network_sent(self):
        """
        测试网络信息采集
        """
        start = time.time()
        current_net_io = psutil.net_io_counters()
        time.sleep(1)

        result = self.collector.get_network_sent_speed(current_net_io.bytes_sent, time.time() - start)
        assert result['name'] == "Network Traffic Sent (KB)"
        assert result['value'] >= 0
        assert result['config'].chart_name == 'Network Traffic (KB)'

    def test_network_recv(self):
        """
        测试网络信息采集
        """
        start = time.time()
        current_net_io = psutil.net_io_counters()
        time.sleep(1)

        result = self.collector.get_network_recv_speed(current_net_io.bytes_recv, time.time() - start)
        assert result['name'] == "Network Traffic Received (KB)"
        assert result['value'] >= 0
        assert result['config'].chart_name == 'Network Traffic (KB)'
