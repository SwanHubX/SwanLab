import numpy as np
import pytest

from swanlab.data.modules.object3d import PointCloud


class TestPointCloud:
    def test_init_validation(self, xyzrgb_points):
        """测试初始化验证"""
        # 正确的初始化
        pc = PointCloud(xyzrgb_points)
        assert pc.points.shape == (4, 6)

        # 错误的形状
        with pytest.raises(ValueError):
            PointCloud(np.zeros((4, 5)))

        # 错误的维度
        with pytest.raises(ValueError):
            PointCloud(np.zeros((4, 4, 3)))

    def test_from_xyz(self, xyz_points):
        """测试从XYZ创建"""
        pc = PointCloud.from_xyz(xyz_points)
        assert pc.points.shape == (4, 6)
        # 验证默认颜色是绿色
        assert np.all(pc.points[:, 3:] == [0, 255, 0])

    def test_from_xyzc(self, xyzc_points):
        """测试从XYZC创建"""
        pc = PointCloud.from_xyzc(xyzc_points)
        assert pc.points.shape == (4, 6)
        # 验证类别被正确映射到颜色
        assert pc.points.shape[1] == 6

    def test_from_xyzrgb(self, xyzrgb_points):
        """测试从XYZRGB创建"""
        pc = PointCloud.from_xyzrgb(xyzrgb_points)
        assert np.array_equal(pc.points, xyzrgb_points)

    def test_metadata(self, xyzrgb_points):
        """测试元数据处理"""
        pc = PointCloud(xyzrgb_points, caption="Test")
        assert pc.caption == "Test"
        assert pc.get_more() == {"caption": "Test"}

    def test_parse(self, xyzrgb_points):
        """测试解析方法"""
        pc = PointCloud(xyzrgb_points)
        filename, buffer = pc.parse()
        assert filename.endswith(".swanlab.pts.json")
        # 默认情况下step为None，所以就生成了None
        assert "stepNone" in filename
        assert buffer.getvalue()  # 确保buffer不为空

    def test_from_swanlab_pts_json_file(self, test_files_dir):
        """测试从文件加载"""
        file_path = test_files_dir / "points.swanlab.pts.json"
        pc = PointCloud.from_swanlab_pts_json_file(file_path)

        # 验证点云数据
        assert pc.points.shape == (4, 6)  # 4个点，每个点6个值

        # 验证具体的点
        expected_points = np.array(
            [
                [0, 0, 0, 255, 0, 0],  # 红色点
                [1, 0, 0, 0, 255, 0],  # 绿色点
                [0, 1, 0, 0, 0, 255],  # 蓝色点
                [0, 0, 1, 255, 255, 255],  # 白色点
            ],
            dtype=np.float32,
        )
        np.testing.assert_array_equal(pc.points, expected_points)

        # 测试错误情况
        with pytest.raises(FileNotFoundError):
            PointCloud.from_swanlab_pts_json_file(test_files_dir / "not_exist.json")

        # 测试带元数据的加载
        pc = PointCloud.from_swanlab_pts_json_file(file_path, step=1, caption="Test Points")
        assert pc.step == 1
        assert pc.caption == "Test Points"
