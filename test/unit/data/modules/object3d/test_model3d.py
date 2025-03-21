from pathlib import Path

import pytest

from swanlab.data.modules.object3d import Model3D


class TestModel3D:
    def test_init_validation(self, test_files_dir):
        """测试初始化验证"""
        glb_path = test_files_dir / "sample.glb"

        # 正确的初始化
        model = Model3D(glb_path)
        assert model.glb_path == glb_path

        # 不存在的文件
        with pytest.raises(FileNotFoundError):
            Model3D(Path("not_exist.glb"))

        # 错误的文件类型
        with pytest.raises(ValueError):
            Model3D(test_files_dir / "wrong.txt")

    def test_from_glb_file(self, test_files_dir):
        """测试从GLB文件创建"""
        glb_path = test_files_dir / "sample.glb"
        model = Model3D.from_glb_file(glb_path, caption="Test")
        assert model.step is None
        assert model.caption == "Test"

    def test_parse(self, test_files_dir):
        """测试解析方法"""
        glb_path = test_files_dir / "sample.glb"
        model = Model3D(glb_path)
        filename, buffer = model.parse()
        assert filename.endswith(".glb")
        # 默认step为None，所以就生成了None
        assert "stepNone" in filename
        assert buffer.getvalue()  # 确保buffer不为空
