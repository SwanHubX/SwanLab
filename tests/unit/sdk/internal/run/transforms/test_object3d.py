"""
@author: cunyue
@file: test_object3d.py
@time: 2026/4/30
@description: Object3D TransformMedia 单元测试
"""

import hashlib
import json

import numpy as np
import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem
from swanlab.sdk.internal.run.transforms.object3d import Object3D

# ---------------------------------- Fixtures ----------------------------------


@pytest.fixture
def xyz_array():
    """(100, 3) XYZ 点云"""
    return np.random.rand(100, 3)


@pytest.fixture
def xyzc_array():
    """(100, 4) XYZC 点云"""
    arr = np.random.rand(100, 4)
    arr[:, 3] = np.random.randint(0, 3, 100)
    return arr


@pytest.fixture
def xyzrgb_array():
    """(100, 6) XYZRGB 点云"""
    arr = np.random.rand(100, 6)
    arr[:, 3:] = np.random.randint(0, 255, (100, 3))
    return arr


@pytest.fixture
def glb_file(tmp_path):
    """生成一个临时 GLB 文件"""
    path = tmp_path / "test.glb"
    path.write_bytes(b"fake glb data for testing")
    return path


@pytest.fixture
def pts_json_file(tmp_path):
    """生成一个临时 .swanlab.pts.json 文件"""
    pts_data = {"version": "0.1", "points": [[0, 0, 0, 0, 255, 0], [1, 1, 1, 255, 0, 0]], "boxes": []}
    path = tmp_path / "test.swanlab.pts.json"
    path.write_text(json.dumps(pts_data))
    return path


# ---------------------------------- 构造测试 ----------------------------------


class TestObject3DInit:
    def test_from_xyz_array(self, xyz_array):
        """从 (N,3) XYZ 数组构造"""
        obj = Object3D(xyz_array)
        assert obj.file_type == "swanlab.pts.json"
        assert obj.caption is None
        assert len(obj.buffer.getvalue()) > 0

    def test_from_xyzc_array(self, xyzc_array):
        """从 (N,4) XYZC 数组构造"""
        obj = Object3D(xyzc_array)
        assert obj.file_type == "swanlab.pts.json"

    def test_from_xyzrgb_array(self, xyzrgb_array):
        """从 (N,6) XYZRGB 数组构造"""
        obj = Object3D(xyzrgb_array)
        assert obj.file_type == "swanlab.pts.json"

    def test_from_dict_points_only(self, xyz_array):
        """从 dict {"points": ...} 构造"""
        obj = Object3D({"points": xyz_array})
        assert obj.file_type == "swanlab.pts.json"

    def test_from_dict_with_boxes(self, xyz_array):
        """从 dict {"points": ..., "boxes": [...]} 构造"""
        boxes = [{"color": [255, 0, 0], "corners": [[0, 0, 0]] * 8, "label": "car", "score": 0.9}]
        obj = Object3D({"points": xyz_array, "boxes": boxes})
        assert obj.file_type == "swanlab.pts.json"
        # 验证 boxes 被序列化到 buffer 中
        data = json.loads(obj.buffer.getvalue())
        assert len(data["boxes"]) == 1
        assert data["boxes"][0]["label"] == "car"

    def test_from_dict_list_points(self):
        """从 dict 的 list 类型 points 构造"""
        points = [[0, 0, 0], [1, 1, 1], [2, 2, 2]]
        obj = Object3D({"points": points})
        assert obj.file_type == "swanlab.pts.json"

    def test_from_glb_file(self, glb_file):
        """从 GLB 文件路径构造"""
        obj = Object3D(str(glb_file))
        assert obj.file_type == "glb"
        assert obj.buffer.getvalue() == b"fake glb data for testing"

    def test_from_glb_path_object(self, glb_file):
        """从 Path 对象构造 GLB"""
        obj = Object3D(glb_file)
        assert obj.file_type == "glb"

    def test_from_pts_json_file(self, pts_json_file):
        """从 .swanlab.pts.json 文件构造"""
        obj = Object3D(str(pts_json_file))
        assert obj.file_type == "swanlab.pts.json"

    def test_caption_stored(self, xyz_array):
        """caption 被正确保存"""
        obj = Object3D(xyz_array, caption="test caption")
        assert obj.caption == "test caption"

    def test_caption_none_by_default(self, xyz_array):
        """caption 默认为 None"""
        obj = Object3D(xyz_array)
        assert obj.caption is None


# ---------------------------------- 构造错误测试 ----------------------------------


class TestObject3DInitErrors:
    def test_unsupported_type_raises(self):
        """传入不支持的类型应抛出 TypeError"""
        with pytest.raises(TypeError, match="Unsupported input type"):
            Object3D(12345)  # type: ignore

    def test_invalid_array_shape_raises(self):
        """不支持的数组形状应抛出 ValueError"""
        with pytest.raises(ValueError, match="Unsupported array shape"):
            Object3D(np.array([1, 2, 3]))

    def test_invalid_channel_count_raises(self):
        """不支持的通道数应抛出 ValueError"""
        with pytest.raises(ValueError, match="Unsupported array shape"):
            Object3D(np.random.rand(100, 5))

    def test_dict_without_points_raises(self):
        """dict 缺少 points 键应抛出 ValueError"""
        with pytest.raises(ValueError, match="must contain 'points' key"):
            Object3D({"boxes": []})

    def test_nonexistent_file_raises(self):
        """不存在的文件应抛出 FileNotFoundError"""
        with pytest.raises(FileNotFoundError, match="File not found"):
            Object3D("nonexistent.glb")

    def test_unsupported_file_type_raises(self, tmp_path):
        """不支持的文件类型应抛出 ValueError"""
        unsupported = tmp_path / "test.xyz"
        unsupported.write_text("data")
        with pytest.raises(ValueError, match="Unsupported file type"):
            Object3D(str(unsupported))

    def test_file_no_extension_raises(self, tmp_path):
        """没有扩展名的文件应抛出 ValueError"""
        no_ext = tmp_path / "noext"
        no_ext.write_text("data")
        with pytest.raises(ValueError, match="no extension"):
            Object3D(str(no_ext))


# ---------------------------------- 套娃加载测试 ----------------------------------


class TestObject3DNesting:
    def test_wrap_copies_buffer(self, xyz_array):
        """套娃加载复用内层 buffer"""
        inner = Object3D(xyz_array, caption="inner")
        outer = Object3D(inner)
        assert outer.buffer is inner.buffer
        assert outer.file_type == inner.file_type
        assert outer.caption == "inner"

    def test_outer_caption_overrides_inner(self, xyz_array):
        """外层 caption 优先级高于内层"""
        inner = Object3D(xyz_array, caption="inner")
        outer = Object3D(inner, caption="outer")
        assert outer.caption == "outer"

    def test_inner_caption_used_when_outer_none(self, xyz_array):
        """外层 caption 为 None 时使用内层 caption"""
        inner = Object3D(xyz_array, caption="inner")
        outer = Object3D(inner, caption=None)
        assert outer.caption == "inner"


# ---------------------------------- column_type 测试 ----------------------------------


class TestObject3DColumnType:
    def test_column_type(self):
        assert Object3D.column_type() == ColumnType.COLUMN_TYPE_OBJECT3D


# ---------------------------------- build_data_record 测试 ----------------------------------


class TestObject3DBuildDataRecord:
    def test_build_data_record_structure(self, xyz_array, tmp_path):
        """build_data_record 返回正确结构的 MediaRecord"""
        obj = Object3D(xyz_array)
        item = obj.transform(step=1, path=tmp_path)
        ts = Timestamp()
        record = Object3D.build_data_record(key="3d", step=1, timestamp=ts, data=[item])
        assert record.key == "3d"
        assert record.step == 1
        assert record.type == ColumnType.COLUMN_TYPE_OBJECT3D
        assert len(record.value.items) == 1

    def test_build_data_record_multiple_items(self, xyz_array, xyzrgb_array, tmp_path):
        """build_data_record 支持多个 MediaItem"""
        i1 = Object3D(xyz_array).transform(step=1, path=tmp_path)
        i2 = Object3D(xyzrgb_array).transform(step=1, path=tmp_path)
        ts = Timestamp()
        record = Object3D.build_data_record(key="k", step=1, timestamp=ts, data=[i1, i2])
        assert len(record.value.items) == 2


# ---------------------------------- transform 测试 ----------------------------------


class TestObject3DTransform:
    def test_transform_returns_media_item(self, xyz_array, tmp_path):
        """transform 返回 MediaItem"""
        item = Object3D(xyz_array).transform(step=1, path=tmp_path)
        assert isinstance(item, MediaItem)

    def test_transform_pts_filename(self, xyz_array, tmp_path):
        """点云 transform 输出 .swanlab.pts.json 文件"""
        item = Object3D(xyz_array).transform(step=1, path=tmp_path)
        assert item.filename.endswith(".swanlab.pts.json")
        assert (tmp_path / item.filename).exists()

    def test_transform_glb_filename(self, glb_file, tmp_path):
        """GLB transform 输出 .glb 文件"""
        item = Object3D(str(glb_file)).transform(step=1, path=tmp_path)
        assert item.filename.endswith(".glb")
        assert (tmp_path / item.filename).exists()

    def test_transform_sha256_correct(self, xyz_array, tmp_path):
        """sha256 与落盘文件内容一致"""
        item = Object3D(xyz_array).transform(step=1, path=tmp_path)
        content = (tmp_path / item.filename).read_bytes()
        assert item.sha256 == hashlib.sha256(content).hexdigest()

    def test_transform_size_correct(self, xyz_array, tmp_path):
        """size 与落盘文件字节数一致"""
        item = Object3D(xyz_array).transform(step=1, path=tmp_path)
        assert item.size == len((tmp_path / item.filename).read_bytes())

    def test_transform_caption_empty_when_none(self, xyz_array, tmp_path):
        """caption 为 None 时返回空字符串"""
        item = Object3D(xyz_array).transform(step=1, path=tmp_path)
        assert item.caption == ""

    def test_transform_caption_preserved(self, xyz_array, tmp_path):
        """caption 有值时正确传递"""
        item = Object3D(xyz_array, caption="hello").transform(step=1, path=tmp_path)
        assert item.caption == "hello"

    def test_pts_json_content_valid(self, xyz_array, tmp_path):
        """点云输出的 JSON 内容合法且包含 version/points/boxes"""
        item = Object3D(xyz_array).transform(step=1, path=tmp_path)
        data = json.loads((tmp_path / item.filename).read_text())
        assert "version" in data
        assert "points" in data
        assert "boxes" in data
        assert len(data["points"]) == 100
        assert len(data["points"][0]) == 6  # XYZRGB


# ---------------------------------- 工厂方法测试 ----------------------------------


class TestObject3DFactoryMethods:
    def test_from_point_data(self, xyz_array):
        """from_point_data 创建 Object3D"""
        obj = Object3D.from_point_data(xyz_array, caption="test")
        assert obj.file_type == "swanlab.pts.json"
        assert obj.caption == "test"

    def test_from_point_data_with_boxes(self, xyz_array):
        """from_point_data 带 boxes"""
        boxes = [{"color": [255, 0, 0], "corners": [[0, 0, 0]] * 8, "label": "car"}]
        obj = Object3D.from_point_data(xyz_array, boxes=boxes, caption="with boxes")
        data = json.loads(obj.buffer.getvalue())
        assert len(data["boxes"]) == 1

    def test_from_point_data_no_boxes(self, xyz_array):
        """from_point_data 不传 boxes 时默认为空"""
        obj = Object3D.from_point_data(xyz_array)
        data = json.loads(obj.buffer.getvalue())
        assert data["boxes"] == []


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
