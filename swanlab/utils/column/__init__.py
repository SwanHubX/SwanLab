from swanlab.sdk.typings.run import ColumnType


def parse_column_type(column: str) -> ColumnType:
    """从前缀中获取指标类型"""
    column_type = column.split(".", 1)[0]
    if column_type == "summary":
        return "SCALAR"
    elif column_type == "config":
        return "CONFIG"
    else:
        return "STABLE"


def to_camel_case(name: str) -> str:
    """将下划线命名转化为驼峰命名"""
    return "".join([w.capitalize() if i > 0 else w for i, w in enumerate(name.split("_"))])
