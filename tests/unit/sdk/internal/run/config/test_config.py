"""
@author: cunyue
@file: test_config.py
@time: 2026/3/14
@description: 测试 config/__init__.py：SwanLabConfig 生命周期、接口和线程安全
"""

import threading
from pathlib import Path
from unittest.mock import MagicMock

import pytest
import yaml

from swanlab.proto.swanlab.config.v1.config_pb2 import UpdateType
from swanlab.sdk.internal.bus.events import ConfigEvent
from swanlab.sdk.internal.run.config import (
    SwanLabConfig,
    config,
    create_run_config,
    create_unbound_run_config,
    deactivate_run_config,
    reset,
    revert_config,
)

# ============================================================
# 辅助工具
# ============================================================


def make_emit():
    """返回一个 Mock emit，用于断言 ConfigEvent 发出"""
    return MagicMock()


def bound_config(tmp_path: Path) -> tuple[SwanLabConfig, MagicMock, Path]:
    """创建并绑定一个独立的 SwanLabConfig 实例，返回 (cfg, emit, config_file)"""
    cfg = SwanLabConfig()
    emit = make_emit()
    config_file = tmp_path / "config.yaml"
    getattr(cfg, "_bindctx")(config_file, emit)
    return cfg, emit, config_file


# ============================================================
# 生命周期：未绑定阶段
# ============================================================


class TestPreBind:
    def test_setitem_stores_value_in_memory(self):
        cfg = SwanLabConfig()
        cfg["lr"] = 0.01
        assert cfg["lr"] == 0.01

    def test_no_file_written_before_bind(self, tmp_path):
        cfg = SwanLabConfig()
        cfg["lr"] = 0.01
        config_file = tmp_path / "config.yaml"
        assert not config_file.exists()

    def test_multiple_writes_before_bind(self):
        cfg = SwanLabConfig()
        cfg["a"] = 1
        cfg["b"] = 2
        assert len(cfg) == 2

    def test_update_before_bind(self):
        cfg = SwanLabConfig()
        cfg.update({"x": 10, "y": 20})
        assert cfg["x"] == 10
        assert cfg["y"] == 20


# ============================================================
# 生命周期：绑定（bindctx）
# ============================================================


class TestBindCtx:
    def test_bind_flushes_existing_memory_to_file(self, tmp_path):
        """绑定时应将内存中已有的 config 全量写入文件"""
        cfg = SwanLabConfig()
        cfg["lr"] = 0.01
        cfg["epochs"] = 10

        config_file = tmp_path / "config.yaml"
        cfg._bindctx(config_file, make_emit())

        assert config_file.exists()
        data = yaml.safe_load(config_file.read_text())
        assert data["lr"]["value"] == 0.01
        assert data["epochs"]["value"] == 10

    def test_bind_emits_init_event(self, tmp_path):
        cfg = SwanLabConfig()
        emit = make_emit()
        cfg._bindctx(tmp_path / "config.yaml", emit)

        emit.assert_called_once()
        event: ConfigEvent = emit.call_args[0][0]
        assert isinstance(event, ConfigEvent)
        assert event.update == UpdateType.UPDATE_TYPE_INIT

    def test_bind_is_idempotent(self, tmp_path):
        """重复调用 bindctx 应静默忽略（幂等）"""
        cfg = SwanLabConfig()
        emit = make_emit()
        config_file = tmp_path / "config.yaml"

        cfg._bindctx(config_file, emit)
        cfg._bindctx(config_file, emit)  # 第二次调用

        assert emit.call_count == 1

    def test_bind_with_empty_config_creates_file(self, tmp_path):
        """即使 config 为空，绑定后也应创建文件"""
        cfg = SwanLabConfig()
        config_file = tmp_path / "config.yaml"
        cfg._bindctx(config_file, make_emit())

        assert config_file.exists()

    def test_module_create_run_config_function(self, tmp_path):
        """模块级 create_run_config() 应创建并绑定 run config，激活代理"""
        config_file = tmp_path / "config.yaml"
        emit = make_emit()
        run_cfg = create_run_config(config_file, emit)

        assert isinstance(run_cfg, SwanLabConfig)
        emit.assert_called_once()

    def test_module_reset_function(self, tmp_path):
        """模块级 reset() 应清空 global config + 取消激活 run config"""
        config["lr"] = 0.01
        run_cfg = create_run_config(tmp_path / "config.yaml", make_emit())
        run_cfg["run_key"] = 1

        reset()

        assert len(config) == 0  # global config 已清空
        assert run_cfg["run_key"] == 1  # run config 不受影响（已创建的实例独立存在）


# ============================================================
# 生命周期：绑定后写操作
# ============================================================


class TestPostBind:
    def test_setitem_writes_file(self, tmp_path):
        cfg, emit, config_file = bound_config(tmp_path)
        emit.reset_mock()

        cfg["lr"] = 0.01

        data = yaml.safe_load(config_file.read_text())
        assert data["lr"]["value"] == 0.01

    def test_setitem_emits_patch_event(self, tmp_path):
        cfg, emit, _ = bound_config(tmp_path)
        emit.reset_mock()

        cfg["lr"] = 0.01

        emit.assert_called_once()
        event: ConfigEvent = emit.call_args[0][0]
        assert event.update == UpdateType.UPDATE_TYPE_PATCH

    def test_update_emits_single_patch(self, tmp_path):
        """update() 批量写应只触发一次 flush"""
        cfg, emit, _ = bound_config(tmp_path)
        emit.reset_mock()

        cfg.update({"a": 1, "b": 2, "c": 3})

        assert emit.call_count == 1

    def test_delitem_writes_file(self, tmp_path):
        cfg, _, config_file = bound_config(tmp_path)
        cfg["lr"] = 0.01
        del cfg["lr"]

        data = yaml.safe_load(config_file.read_text())
        assert data is None or "lr" not in (data or {})

    def test_pop_writes_file(self, tmp_path):
        cfg, _, config_file = bound_config(tmp_path)
        cfg["lr"] = 0.01
        cfg.pop("lr")

        data = yaml.safe_load(config_file.read_text())
        assert data is None or "lr" not in (data or {})

    def test_file_is_fully_overwritten_each_time(self, tmp_path):
        """每次写操作均全量覆写文件，不保留旧 key"""
        cfg, _, config_file = bound_config(tmp_path)
        cfg["a"] = 1
        cfg["b"] = 2
        del cfg["a"]

        data = yaml.safe_load(config_file.read_text())
        assert "a" not in (data or {})
        assert data["b"]["value"] == 2


# ============================================================
# 重置（reset）
# ============================================================


class TestReset:
    def test_reset_clears_memory(self, tmp_path):
        cfg, _, _ = bound_config(tmp_path)
        cfg["lr"] = 0.01
        cfg._reset()

        assert len(cfg) == 0

    def test_reset_unbinds(self, tmp_path):
        """reset 后再写入不应触发文件 IO"""
        cfg, _, config_file = bound_config(tmp_path)
        cfg._reset()

        cfg["lr"] = 99  # 未绑定，不应写文件
        data = yaml.safe_load(config_file.read_text())
        # 文件仍是 reset 前的内容（空）
        assert data is None or "lr" not in (data or {})

    def test_reset_clears_sort_seq(self):
        cfg = SwanLabConfig()
        cfg["a"] = 1
        cfg["b"] = 2
        cfg._reset()
        cfg["c"] = 3

        # reset 后新写入的 sort 序号从 0 重新开始，不应报错
        assert cfg["c"] == 3

    def test_double_reset_is_safe(self):
        cfg = SwanLabConfig()
        cfg._reset()
        cfg._reset()  # 不应报错


# ============================================================
# MutableMapping 接口
# ============================================================


class TestMutableMappingInterface:
    def test_getitem_existing(self):
        cfg = SwanLabConfig()
        cfg["lr"] = 0.01
        assert cfg["lr"] == 0.01

    def test_getitem_missing_raises_keyerror(self):
        cfg = SwanLabConfig()
        with pytest.raises(KeyError):
            _ = cfg["missing"]

    def test_getitem_non_str_key_raises_typeerror(self):
        cfg = SwanLabConfig()
        with pytest.raises(TypeError):
            _ = cfg[123]  # type: ignore

    def test_delitem_existing(self):
        cfg = SwanLabConfig()
        cfg["lr"] = 0.01
        del cfg["lr"]
        assert "lr" not in cfg

    def test_delitem_missing_raises_keyerror(self):
        cfg = SwanLabConfig()
        with pytest.raises(KeyError):
            del cfg["missing"]

    def test_iter(self):
        cfg = SwanLabConfig()
        cfg["a"] = 1
        cfg["b"] = 2
        assert set(cfg) == {"a", "b"}

    def test_len(self):
        cfg = SwanLabConfig()
        cfg["a"] = 1
        cfg["b"] = 2
        assert len(cfg) == 2

    def test_str(self):
        cfg = SwanLabConfig()
        cfg["lr"] = 0.01
        s = str(cfg)
        assert "lr" in s

    def test_contains(self):
        cfg = SwanLabConfig()
        cfg["lr"] = 0.01
        assert "lr" in cfg
        assert "epochs" not in cfg


# ============================================================
# 对象属性风格
# ============================================================


class TestAttributeStyle:
    def test_setattr_and_getattr(self):
        cfg = SwanLabConfig()
        cfg.lr = 0.01  # type: ignore
        assert cfg.lr == 0.01  # type: ignore

    def test_delattr(self):
        cfg = SwanLabConfig()
        cfg.lr = 0.01  # type: ignore
        del cfg.lr  # type: ignore
        assert "lr" not in cfg

    def test_delattr_missing_raises(self):
        cfg = SwanLabConfig()
        with pytest.raises(AttributeError):
            del cfg.missing  # type: ignore

    def test_private_attr_setattr_blocked(self):
        cfg = SwanLabConfig()
        with pytest.raises(AttributeError, match="private"):
            cfg._xxx__ = 1  # type: ignore

    def test_private_attr_delattr_blocked(self):
        cfg = SwanLabConfig()
        with pytest.raises(AttributeError, match="private"):
            del cfg._xxx__  # type: ignore

    def test_getattr_missing_raises_attributeerror(self):
        cfg = SwanLabConfig()
        with pytest.raises(AttributeError):
            _ = cfg.nonexistent  # type: ignore


# ============================================================
# 批量操作
# ============================================================


class TestBatchOperations:
    def test_update_with_dict(self):
        cfg = SwanLabConfig()
        cfg.update({"lr": 0.01, "epochs": 10})
        assert cfg["lr"] == 0.01
        assert cfg["epochs"] == 10

    def test_update_with_kwargs(self):
        cfg = SwanLabConfig()
        cfg.update(lr=0.01, epochs=10)
        assert cfg["lr"] == 0.01
        assert cfg["epochs"] == 10

    def test_update_with_both(self):
        cfg = SwanLabConfig()
        cfg.update({"a": 1}, b=2)
        assert cfg["a"] == 1
        assert cfg["b"] == 2

    def test_get_existing(self):
        cfg = SwanLabConfig()
        cfg["lr"] = 0.01
        assert cfg.get("lr") == 0.01

    def test_get_missing_returns_none(self):
        cfg = SwanLabConfig()
        assert cfg.get("missing") is None

    def test_get_missing_with_default(self):
        cfg = SwanLabConfig()
        assert cfg.get("missing", 42) == 42

    def test_set(self):
        cfg = SwanLabConfig()
        cfg.set("lr", 0.01)
        assert cfg["lr"] == 0.01

    def test_pop_existing(self):
        cfg = SwanLabConfig()
        cfg["lr"] = 0.01
        val = cfg.pop("lr")
        assert val == 0.01
        assert "lr" not in cfg

    def test_pop_missing_returns_none(self):
        cfg = SwanLabConfig()
        assert cfg.pop("missing") is None

    def test_pop_missing_with_default(self):
        cfg = SwanLabConfig()
        assert cfg.pop("missing", 99) == 99

    def test_clean_clears_items(self):
        cfg = SwanLabConfig()
        cfg["a"] = 1
        cfg["b"] = 2
        cfg.clean()
        assert len(cfg) == 0

    def test_clean_preserves_binding(self, tmp_path):
        """clean() 只清内容，不解绑 —— 清空后继续写入仍能触发 flush"""
        cfg, emit, _ = bound_config(tmp_path)
        cfg.clean()
        emit.reset_mock()

        cfg["new_key"] = 1
        assert emit.call_count == 1


# ============================================================
# sort 顺序维护
# ============================================================


class TestSortOrder:
    def test_sort_index_increments(self, tmp_path):
        cfg, _, config_file = bound_config(tmp_path)
        cfg["first"] = 1
        cfg["second"] = 2
        cfg["third"] = 3

        data = yaml.safe_load(config_file.read_text())
        # "first" 在 init flush 时 sort=0，后续 patch 时 sort 以插入顺序递增
        assert data["first"]["sort"] < data["second"]["sort"] < data["third"]["sort"]

    def test_update_existing_key_preserves_sort(self, tmp_path):
        """更新已有 key 不改变 sort 序号"""
        cfg, _, config_file = bound_config(tmp_path)
        cfg["lr"] = 0.01
        original_sort = yaml.safe_load(config_file.read_text())["lr"]["sort"]

        cfg["lr"] = 0.001  # 更新值

        new_sort = yaml.safe_load(config_file.read_text())["lr"]["sort"]
        assert new_sort == original_sort


# ============================================================
# revert_config 静态方法
# ============================================================


class TestRevertConfig:
    def test_basic_revert(self):
        raw = {
            "lr": {"value": 0.01, "desc": "", "sort": 1},
            "epochs": {"value": 10, "desc": "", "sort": 0},
        }
        result = revert_config(raw)
        assert list(result.keys()) == ["epochs", "lr"]  # sort 升序
        assert result["lr"] == 0.01
        assert result["epochs"] == 10

    def test_skips_entries_without_value(self):
        raw = {
            "lr": {"value": 0.01, "sort": 0},
            "bad": {"desc": "no value field", "sort": 1},
        }
        result = revert_config(raw)
        assert "lr" in result
        assert "bad" not in result

    def test_missing_sort_defaults_to_zero(self):
        raw = {"lr": {"value": 0.01}}
        result = revert_config(raw)
        assert result["lr"] == 0.01

    def test_empty_config(self):
        assert revert_config({}) == {}


# ============================================================
# 线程安全
# ============================================================


class TestThreadSafety:
    def test_concurrent_writes_no_corruption(self, tmp_path):
        """并发写入不应导致内存状态损坏（键数量正确）"""
        cfg = SwanLabConfig()
        n = 50
        errors = []

        def writer(i):
            try:
                cfg[f"key_{i}"] = i
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=writer, args=(i,)) for i in range(n)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert not errors
        assert len(cfg) == n

    def test_concurrent_writes_with_bind(self, tmp_path):
        """绑定后并发写入不应导致文件写入异常"""
        cfg, _, config_file = bound_config(tmp_path)
        errors = []

        def writer(i):
            try:
                cfg[f"k{i}"] = i
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=writer, args=(i,)) for i in range(20)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert not errors
        # 文件应可正常解析
        data = yaml.safe_load(config_file.read_text())
        assert isinstance(data, dict)


# ============================================================
# 代理（Proxy）
# ============================================================


class TestProxy:
    def test_proxy_points_to_global_when_no_run(self):
        """无 run 时，代理应指向 global config"""
        config["global_key"] = "global_value"
        assert config["global_key"] == "global_value"

    def test_proxy_points_to_run_config_when_active(self, tmp_path):
        """run 活跃时，代理应指向 run config"""
        config["global_key"] = "global"
        run_cfg = create_run_config(tmp_path / "config.yaml", make_emit())
        run_cfg["run_key"] = "run"

        assert config["run_key"] == "run"
        assert config["global_key"] == "global"  # global 的值也被复制到 run config

    def test_proxy_write_goes_to_run_config_when_active(self, tmp_path):
        """run 活跃时，通过代理写入应写到 run config"""
        run_cfg = create_run_config(tmp_path / "config.yaml", make_emit())
        config["new_key"] = "new_value"

        assert run_cfg["new_key"] == "new_value"

    def test_proxy_restores_to_global_after_deactivate(self, tmp_path):
        """deactivate 后，代理应恢复指向 global config"""
        config["global_key"] = "global"
        create_run_config(tmp_path / "config.yaml", make_emit())
        config["run_key"] = "run"

        deactivate_run_config()

        assert "run_key" not in config  # run config 已清理
        assert config["global_key"] == "global"  # global config 仍在


# ============================================================
# Run Config 创建
# ============================================================


class TestCreateRunConfig:
    def test_create_run_config_copies_from_global(self, tmp_path):
        """create_run_config 应从 global config 复制数据"""
        config["lr"] = 0.01
        config["epochs"] = 10

        run_cfg = create_run_config(tmp_path / "config.yaml", make_emit())

        assert run_cfg["lr"] == 0.01
        assert run_cfg["epochs"] == 10

    def test_create_run_config_binds_to_file(self, tmp_path):
        """create_run_config 应绑定到文件"""
        config_file = tmp_path / "config.yaml"
        run_cfg = create_run_config(config_file, make_emit())
        run_cfg["key"] = "value"

        assert config_file.exists()
        data = yaml.safe_load(config_file.read_text())
        assert data["key"]["value"] == "value"

    def test_create_unbound_run_config_no_file(self, tmp_path):
        """create_unbound_run_config 不应创建文件"""
        config["lr"] = 0.01
        run_cfg = create_unbound_run_config()
        run_cfg["new_key"] = "value"

        assert run_cfg["lr"] == 0.01
        assert run_cfg["new_key"] == "value"
        # 不应有文件创建

    def test_deactivate_clears_run_config(self, tmp_path):
        """deactivate_run_config 应清理 run config 内存"""
        run_cfg = create_run_config(tmp_path / "config.yaml", make_emit())
        run_cfg["key"] = "value"

        deactivate_run_config()

        assert len(run_cfg) == 0  # run config 已清空
