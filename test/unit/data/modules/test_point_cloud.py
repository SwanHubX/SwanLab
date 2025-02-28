import tempfile
import os
import json
import numpy as np
import pytest
from swanlab.data.modules import PointCloud
from swanlab.data.modules.point_cloud import PointsData


def test_point_cloud_xyz():
    """测试xyz格式点云"""
    # 生成随机xyz点云数据
    points = np.random.rand(100, 3) * 100  # 100个点
    pc = PointCloud(points)

    # 测试基本属性
    assert pc.points.format == 'xyz'
    assert pc.points.take().shape == (100, 6)  # 转换后变为xyzrgb格式

    # 测试parse结果
    filename, buffer = pc.parse()
    assert isinstance(filename, str)
    assert filename.startswith("pointscloud-step")
    assert filename.endswith(".json")
    assert buffer is not None

    # 测试caption为None的情况
    assert pc.get_more() is None

    # 测试分组名和图表类型
    assert pc.get_section() == "PointsCloud"
    assert pc.get_chart() == pc.Chart.ChartItem


def test_point_cloud_xyzc():
    """测试xyzc格式点云"""
    # 生成带类别的点云数据
    xyz = np.random.rand(100, 3) * 100
    categories = np.random.randint(0, 15, (100, 1))
    points = np.hstack([xyz, categories])

    pc = PointCloud(points)
    assert pc.points.format == 'xyzc'
    assert pc.points.take().shape == (100, 6)  # 转换后变为xyzrgb格式

    # 测试颜色映射是否正确
    colors = pc.points.take()[:, 3:]
    assert np.all((colors >= 0) & (colors <= 255))


def test_point_cloud_xyzrgb():
    """测试xyzrgb格式点云"""
    # 生成带RGB的点云数据
    xyz = np.random.rand(100, 3) * 100
    rgb = np.random.randint(0, 255, (100, 3))
    points = np.hstack([xyz, rgb])

    pc = PointCloud(points)
    assert pc.points.format == 'xyzrgb'
    assert pc.points.take().shape == (100, 6)

    # RGB值应保持不变
    np.testing.assert_array_equal(pc.points.take()[:, 3:], rgb)


def test_point_cloud_with_caption():
    """测试带caption的点云"""
    points = np.random.rand(100, 3) * 100
    caption = "Test Point Cloud"
    pc = PointCloud(points, caption=caption)

    more_info = pc.get_more()
    assert more_info is not None
    assert more_info["caption"] == caption


def test_point_cloud_validation():
    """测试点云数据验证"""
    # 测试None输入
    with pytest.raises(ValueError, match="Points array cannot be None"):
        PointCloud(None)

    # 测试空数组
    with pytest.raises(ValueError, match="Points array cannot be empty"):
        PointCloud(np.array([]))

    # 测试维度错误
    with pytest.raises(ValueError, match="Points must be a 2D array"):
        PointCloud(np.random.rand(100, 3, 3))

    # 测试特征数量错误
    with pytest.raises(ValueError, match="Point cloud data must be in one of these formats"):
        PointCloud(np.random.rand(100, 5))


def test_point_cloud_category_validation():
    """测试类别验证"""
    # 生成非整数类别
    xyz = np.random.rand(100, 3)
    categories = np.random.rand(100, 1) * 20  # 非整数类别
    points = np.hstack([xyz, categories])

    with pytest.raises(ValueError, match="Categories must be integers"):
        PointCloud(points)

    # 类别超出范围
    categories = np.random.randint(16, 20, (100, 1))  # 超出0-15范围
    points = np.hstack([xyz, categories])

    with pytest.raises(ValueError, match="Categories must be in range"):
        PointCloud(points)


def test_point_cloud_rgb_validation():
    """测试RGB验证"""
    # 生成非整数RGB值
    xyz = np.random.rand(100, 3)
    rgb = np.random.rand(100, 3) * 300  # 非整数RGB
    points = np.hstack([xyz, rgb])

    with pytest.raises(ValueError, match="RGB values must be integers"):
        PointCloud(points)

    # RGB值超出范围
    rgb = np.random.randint(256, 300, (100, 3))  # 超出0-255范围
    points = np.hstack([xyz, rgb])

    with pytest.raises(ValueError, match="RGB values must be in range"):
        PointCloud(points)


def test_point_cloud_factory_methods():
    """测试工厂方法"""
    xyz = np.random.rand(100, 3)
    xyz_pc = PointsData.from_xyz(xyz)
    assert xyz_pc.format == 'xyz'

    xyzc = np.hstack([xyz, np.zeros((100, 1))])
    xyzc_pc = PointsData.from_xyzc(xyzc)
    assert xyzc_pc.format == 'xyzc'

    xyzrgb = np.hstack([xyz, np.zeros((100, 3))])
    xyzrgb_pc = PointsData.from_xyzrgb(xyzrgb)
    assert xyzrgb_pc.format == 'xyzrgb'

    # 测试错误的输入
    with pytest.raises(ValueError, match="XYZ data must have 3 features"):
        PointsData.from_xyz(xyzc)

    with pytest.raises(ValueError, match="XYZC data must have 4 features"):
        PointsData.from_xyzc(xyz)

    with pytest.raises(ValueError, match="XYZRGB data must have 6 features"):
        PointsData.from_xyzrgb(xyzc)


def test_point_cloud_data_types():
    """测试不同数据类型"""
    dtypes = [np.float32, np.float64, np.int32, np.int64]
    for dtype in dtypes:
        points = np.random.rand(100, 3).astype(dtype)
        pc = PointCloud(points)
        assert pc.points.take().dtype == dtype


@pytest.mark.parametrize("n_points", [1, 10, 1000])
def test_point_cloud_different_sizes(n_points):
    """测试不同大小的点云"""
    points = np.random.rand(n_points, 3)
    pc = PointCloud(points)
    assert pc.points.take().shape == (n_points, 6)


def test_point_cloud_from_json_file():
    """测试从JSON文件加载点云数据的正常情况"""
    # 创建临时文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        # 测试xyz格式
        points_xyz = np.random.rand(100, 3) * 100
        json.dump(points_xyz.tolist(), f)
        f.flush()

        # 测试加载
        pc = PointCloud.from_json_file(f.name, caption="Test XYZ")
        assert pc.points.format == 'xyz'
        np.testing.assert_array_equal(pc.points.take()[:, :3], points_xyz)
        assert pc.caption == "Test XYZ"

    # 清理临时文件
    os.unlink(f.name)


def test_point_cloud_from_json_file_xyzc():
    """测试从JSON文件加载xyzc格式点云数据"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        xyz = np.random.rand(100, 3) * 100
        categories = np.random.randint(0, 15, (100, 1))
        points_xyzc = np.hstack([xyz, categories])
        json.dump(points_xyzc.tolist(), f)
        f.flush()

        pc = PointCloud.from_json_file(f.name)
        assert pc.points.format == 'xyzc'
        np.testing.assert_array_equal(pc.points.take()[:, :3], xyz)

    os.unlink(f.name)


def test_point_cloud_from_json_file_xyzrgb():
    """测试从JSON文件加载xyzrgb格式点云数据"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        xyz = np.random.rand(100, 3) * 100
        rgb = np.random.randint(0, 255, (100, 3))
        points_xyzrgb = np.hstack([xyz, rgb])
        json.dump(points_xyzrgb.tolist(), f)
        f.flush()

        pc = PointCloud.from_json_file(f.name)
        assert pc.points.format == 'xyzrgb'
        np.testing.assert_array_equal(pc.points.take(), points_xyzrgb)

    os.unlink(f.name)


def test_point_cloud_from_json_file_errors():
    """测试从JSON文件加载时的错误情况"""
    # 测试文件不存在
    with pytest.raises(FileNotFoundError):
        PointCloud.from_json_file("nonexistent_file.json")

    # 测试无效的JSON格式
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        f.write("invalid json content")
        f.flush()

        with pytest.raises(ValueError, match="Invalid JSON format"):
            PointCloud.from_json_file(f.name)

    os.unlink(f.name)

    # 测试非法的点云数据格式
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        # 写入非嵌套列表
        json.dump([1, 2, 3], f)
        f.flush()

        with pytest.raises(ValueError, match="JSON file must contain a list of point lists"):
            PointCloud.from_json_file(f.name)

    os.unlink(f.name)

    # 测试错误的维度
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        # 写入错误维度的数据
        points = np.random.rand(100, 5)  # 5维数据是不支持的
        json.dump(points.tolist(), f)
        f.flush()

        with pytest.raises(ValueError, match="Point cloud data must be in one of these formats"):
            PointCloud.from_json_file(f.name)

    os.unlink(f.name)


def test_point_cloud_from_json_file_validation():
    """测试从JSON文件加载时的数据验证"""
    # 测试非法的类别值
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        xyz = np.random.rand(100, 3)
        categories = np.random.randint(16, 20, (100, 1))  # 超出范围的类别
        points = np.hstack([xyz, categories])
        json.dump(points.tolist(), f)
        f.flush()

        with pytest.raises(ValueError, match="Categories must be in range"):
            PointCloud.from_json_file(f.name)

    os.unlink(f.name)

    # 测试非法的RGB值
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        xyz = np.random.rand(100, 3)
        rgb = np.random.randint(256, 300, (100, 3))  # 超出范围的RGB值
        points = np.hstack([xyz, rgb])
        json.dump(points.tolist(), f)
        f.flush()

        with pytest.raises(ValueError, match="RGB values must be in range"):
            PointCloud.from_json_file(f.name)

    os.unlink(f.name)
