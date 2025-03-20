import json

import numpy as np
import pytest


@pytest.fixture
def xyz_points():
    """生成示例XYZ点云数据"""
    return np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ],
        dtype=np.float32,
    )


@pytest.fixture
def xyzc_points():
    """生成示例XYZC点云数据"""
    points = np.zeros((4, 4), dtype=np.float32)
    points[:, :3] = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    points[:, 3] = [0, 1, 2, 3]  # 类别标签
    return points


@pytest.fixture
def xyzrgb_points():
    """生成示例XYZRGB点云数据"""
    points = np.zeros((4, 6), dtype=np.float32)
    points[:, :3] = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]
    )
    points[:, 3:] = np.array(
        [
            [255, 0, 0],
            [0, 255, 0],
            [0, 0, 255],
            [255, 255, 255],
        ]
    )
    return points


@pytest.fixture(scope="session")
def test_files_dir(tmp_path_factory):
    """创建会话级别的测试文件"""
    test_dir = tmp_path_factory.mktemp("test_files")

    # 创建示例GLB文件
    glb_file = test_dir / "sample.glb"
    glb_file.write_bytes(b'dummy glb content')

    # 创建错误类型文件
    txt_file = test_dir / "wrong.txt"
    txt_file.write_text('dummy text content')

    # 创建示例点云文件
    points_data = [
        [0, 0, 0, 255, 0, 0],  # 红色点
        [1, 0, 0, 0, 255, 0],  # 绿色点
        [0, 1, 0, 0, 0, 255],  # 蓝色点
        [0, 0, 1, 255, 255, 255],  # 白色点
    ]
    data = {"points": points_data}
    pts_file = test_dir / "points.swanlab.pts.json"
    pts_file.write_text(json.dumps(data))

    return test_dir
