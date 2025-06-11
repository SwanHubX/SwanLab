import json

import pytest

from swanlab.data.modules.custom_charts.table import Table


def test_table():
    """测试表格功能"""
    # 创建表格对象
    table = Table()
    assert isinstance(table, Table)

    # 添加数据
    headers = ["姓名", "年龄", "城市"]
    rows = [["张三", 25, "北京"], ["李四", 30, "上海"]]
    table.add(headers, rows)

    # 不允许设置标题
    with pytest.raises(NotImplementedError):
        table.set_global_opts({"title": "测试表格"})

    # 测试格式转换
    formatted_data = table.get_table_format()
    assert isinstance(formatted_data, dict)
    assert "_swanLab" in formatted_data
    assert formatted_data["_swanLab"] == "table"
    assert "colDefs" in formatted_data
    assert "rowData" in formatted_data

    # 使用标准json序列化测试数据结构
    json_str = json.dumps(formatted_data)
    assert isinstance(json_str, str)

    # 可解析为字典
    table_data = json.loads(json_str)
    assert isinstance(table_data, dict)
    assert table_data["_swanLab"] == "table"

    # 验证列定义
    assert len(table_data["colDefs"]) == 3
    assert table_data["colDefs"][0]["field"] == "姓名"

    # 验证行数据
    assert len(table_data["rowData"]) == 2
    assert table_data["rowData"][0]["姓名"] == "张三"
    assert table_data["rowData"][0]["年龄"] == 25
