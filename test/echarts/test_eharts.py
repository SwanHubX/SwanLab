import swanlab

class TestBarChart:
    """测试柱状图基础功能"""

    def test_create_bar_chart(self):
        """测试创建柱状图"""
        chart = (
            swanlab.echarts.Bar()
            .add_xaxis(["A", "B", "C"])
            .add_yaxis("数据", [1, 2, 3])
        )
        assert chart is not None


class TestChartLogging:
    """测试图表打印功能"""

    def setup_method(self):
        """测试前初始化"""
        swanlab.init(experiment_name="test_charts", mode="disabled")

    def test_log_chart(self):
        """测试记录图表"""
        chart = (
            swanlab.echarts.Bar()
            .add_xaxis(["X", "Y"])
            .add_yaxis("数值", [10, 20])
        )
        swanlab.log({"chart": chart})