import numpy as np
import pytest

from swanlab.data.modules.object3d import Model3D, Object3D, PointCloud


class TestObject3D:
    def test_from_xyz_array(self, xyz_points):
        """测试从XYZ数组创建"""
        obj = Object3D(xyz_points)
        assert isinstance(obj, PointCloud)
        assert obj.points.shape == (4, 6)

    def test_from_xyzc_array(self, xyzc_points):
        """测试从XYZC数组创建"""
        obj = Object3D(xyzc_points)
        assert isinstance(obj, PointCloud)
        assert obj.points.shape == (4, 6)

    def test_from_xyzrgb_array(self, xyzrgb_points):
        """测试从XYZRGB数组创建"""
        obj = Object3D(xyzrgb_points)
        assert isinstance(obj, PointCloud)
        assert np.array_equal(obj.points, xyzrgb_points)

    def test_from_glb_file(self, test_files_dir):
        """测试从GLB文件创建"""
        glb_path = test_files_dir / "sample.glb"
        obj = Object3D(glb_path)
        assert isinstance(obj, Model3D)

    def test_from_swanlab_pts_json(self, test_files_dir):
        """测试从SwanLab点云文件创建"""
        pts_path = test_files_dir / "points.swanlab.pts.json"
        obj = Object3D(pts_path)
        assert isinstance(obj, PointCloud)

    def test_unsupported_format(self, test_files_dir):
        """测试不支持的格式"""
        with pytest.raises(ValueError):
            Object3D(np.zeros((4, 5)))  # 错误的点云形状

        # 测试不支持的文件类型
        wrong_file = test_files_dir / "wrong.txt"
        with pytest.raises(ValueError, match="Unsupported file type"):
            Object3D(wrong_file)

    def test_metadata_handling(self, xyz_points):
        """测试元数据处理"""
        obj = Object3D(xyz_points, caption="Test")
        assert obj.step is None
        assert obj.caption == "Test"

    def test_nested_mediatype(self, xyz_points):
        """测试Object3D类支持嵌套MediaType输入（套娃）"""
        # 创建基础Object3D实例
        base_obj = Object3D(xyz_points, caption="original")
        
        # 测试嵌套输入 - 直接返回原实例
        nested_obj = Object3D(base_obj)
        assert nested_obj is base_obj  # 应该返回同一个实例
        
        # 测试可以覆盖caption
        nested_with_new_caption = Object3D(base_obj, caption="new caption")
        # 对于PointCloud等类型，需要检查是否创建了新实例或更新了caption
        assert nested_with_new_caption.caption == "new caption"
        
        # 测试从Molecule创建的嵌套
        try:
            from rdkit.Chem import MolFromSmiles
            mol = MolFromSmiles("CCO")
            if mol:
                base_molecule = Object3D(mol, caption="molecule original")
                nested_molecule = Object3D(base_molecule)
                # Object3D应该返回Molecule实例
                assert nested_molecule.caption == "molecule original"
                
                nested_molecule_new_caption = Object3D(base_molecule, caption="molecule new")
                assert nested_molecule_new_caption.caption == "molecule new"
        except ImportError:
            # RDKit不可用时跳过Molecule测试
            pass
