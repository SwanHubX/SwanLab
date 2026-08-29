"""DefinitionResolver 单元测试。

覆盖：glob 路由（glob 形式合法性由 Run.define_metric 在 API 层校验，
见 test_run_define_metric.py；事件经 ``is_glob`` 显式携带分类结果）、
metric class 适配（MEDIA 回退 _step、SCALAR 自引用放行）、解析优先级
（exact > 最长前缀 glob > automatic 钉住）、X 值去重、
merge/replace 语义（含 step_sync 默认值 = is_custom_x(x_axis)）。
"""

import pytest

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.sdk.internal.bus.events import MetricDefineEvent
from swanlab.sdk.internal.run.components.consumer.resolver import DefinitionResolver


def _define(resolver: DefinitionResolver, key: str, is_glob: bool = False) -> None:
    """对给定 key 触发一次 handle_define（其余字段默认未提供）。"""
    resolver.handle_define(MetricDefineEvent(key=key, is_glob=is_glob))


class TestGlobRouting:
    """事件 is_glob 的 exact/glob 分区路由。"""

    def test_exact_key_accepted(self):
        """无 '*' 的 key 作为 exact rule 注册。"""
        r = DefinitionResolver()
        _define(r, "train/loss")
        assert "train/loss" in r._exact_rules
        assert len(r._glob_rules) == 0

    @pytest.mark.parametrize("pattern", ["train/*", "*"], ids=lambda v: repr(v))
    def test_glob_patterns_route_to_glob_rules(self, pattern):
        """末尾单 '*'（含单独的 '*'）作为 glob rule 注册。"""
        r = DefinitionResolver()
        _define(r, pattern, is_glob=True)
        assert pattern in r._glob_rules
        assert len(r._exact_rules) == 0


class TestSelfReferenceXAxis:
    def test_scalar_self_reference_x_axis_is_kept(self):
        """SCALAR 允许 x_axis == key（图像为直线 y=x），不再降级为 _step。"""
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("epoch", x_axis="epoch"))
        concrete = r.resolve_concrete("epoch", "SCALAR")
        assert concrete.effective.x_axis == "epoch"
        column = r.materialize_column("epoch", "SCALAR", ColumnType.COLUMN_TYPE_SCALAR)
        assert column is not None
        assert column.x_axis == "epoch"

    def test_media_x_axis_still_falls_back_to_step(self):
        """MEDIA 一律回退 _step（既有行为不变）。"""
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("train/img", x_axis="epoch"))
        concrete = r.resolve_concrete("train/img", "MEDIA")
        assert concrete.effective.x_axis == "_step"


class TestWandbAlignedResolution:
    def test_glob_update_does_not_retroactively_change_materialized_concrete(self):
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("train/*", is_glob=True, section_name="A"))
        first = r.resolve_concrete("train/loss", "SCALAR")
        r.materialize_column("train/loss", "SCALAR", ColumnType.COLUMN_TYPE_SCALAR)

        r.handle_define(MetricDefineEvent("train/*", is_glob=True, section_name="B"))

        assert r.resolve_concrete("train/loss", "SCALAR") is first
        assert first.effective.section_name == "A"
        assert r.resolve_concrete("train/new_loss", "SCALAR").effective.section_name == "B"

    def test_exact_define_does_not_retroactively_change_glob_concrete(self):
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("train/*", is_glob=True, section_name="Train", hidden=True))
        r.resolve_concrete("train/loss", "SCALAR")
        r.materialize_column("train/loss", "SCALAR", ColumnType.COLUMN_TYPE_SCALAR)

        columns = r.handle_define(MetricDefineEvent("train/loss", x_axis="epoch"))

        # 图表定义 first-writer-wins：exact define 不再对已物化 concrete 产出 ColumnRecord
        assert columns is None
        # rule 仍登记，供之后首次出现的 key 使用
        assert "train/loss" in r._exact_rules
        # 已物化的 glob concrete 不被回溯（保持 glob 快照）
        effective = r.resolve_concrete("train/loss", "SCALAR").effective
        assert effective.x_axis == "_step"
        assert effective.section_name == "Train"
        assert effective.hidden is True

    def test_multiple_defines_before_log_last_define_wins(self):
        """同一 key 在 log 前多次 define，最后一次定义生效（merge 累积在 rule 上）。"""
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("train/loss", x_axis="epoch"))
        r.handle_define(MetricDefineEvent("train/loss", x_axis="iter", section_name="Training"))

        concrete = r.resolve_concrete("train/loss", "SCALAR")
        assert concrete.effective.x_axis == "iter"
        assert concrete.effective.section_name == "Training"

        column = r.materialize_column("train/loss", "SCALAR", ColumnType.COLUMN_TYPE_SCALAR)
        assert column is not None
        assert column.x_axis == "iter"
        assert column.section_name == "Training"

    def test_hidden_three_state_merge_semantics(self):
        """hidden 三态：None=保留旧值，True=隐藏，False=显式解除（merge 模式下同样生效）。"""
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("train/loss", hidden=True))
        # 未提供（None）→ 保留 True
        r.handle_define(MetricDefineEvent("train/loss", section_name="Train"))
        assert r._exact_rules["train/loss"].hidden is True
        # 显式 False → merge 下解除隐藏
        r.handle_define(MetricDefineEvent("train/loss", hidden=False))
        assert r._exact_rules["train/loss"].hidden is False

    def test_hidden_replace_resets_to_default(self):
        """overwrite=True（replace）下未提供 hidden → 重置为 False。"""
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("train/loss", hidden=True))
        r.handle_define(MetricDefineEvent("train/loss", overwrite=True))
        assert r._exact_rules["train/loss"].hidden is False

    def test_subsequent_log_after_define_does_not_upgrade_automatic_concrete(self):
        """log → define → 再 log：concrete 保持首次 automatic 快照，列始终由 core 自动建。"""
        r = DefinitionResolver()
        # ① 首次 log（无 rule → automatic 钉住，不产出列）
        first = r.resolve_concrete("train/loss", "SCALAR")
        assert first.source == "automatic"
        col1 = r.materialize_column("train/loss", "SCALAR", ColumnType.COLUMN_TYPE_SCALAR)
        assert col1 is None  # automatic 不产出 ColumnRecord，列由 core 收到数据后自动创建

        # ② define 登记 exact rule（x_axis=epoch），不产出 ColumnRecord
        columns = r.handle_define(MetricDefineEvent("train/loss", x_axis="epoch"))
        assert columns is None
        assert "train/loss" in r._exact_rules  # rule 仍登记

        # ③ 再 log：resolve_concrete 返回冻结的首次 concrete，不升级为 exact
        second = r.resolve_concrete("train/loss", "SCALAR")
        assert second is first  # 同一对象
        assert second.source == "automatic"  # 仍是 automatic，未被 exact 认领
        assert second.effective.x_axis == "_step"  # x_axis 未变

        # ④ materialize_column 仍不产出（automatic 来源永不物化）
        col2 = r.materialize_column("train/loss", "SCALAR", ColumnType.COLUMN_TYPE_SCALAR)
        assert col2 is None


class TestXValueDedup:
    """try_accept_x_value：custom X 值 first-writer-wins 去重。"""

    def test_new_value_accepted(self):
        """首次出现的值与其后的不同值均被接受。"""
        r = DefinitionResolver()
        assert r.try_accept_x_value("acc", 5.0) is True
        assert r.try_accept_x_value("acc", 6.0) is True

    def test_same_value_rejected(self):
        r = DefinitionResolver()
        r.try_accept_x_value("acc", 5.0)
        assert r.try_accept_x_value("acc", 5.0) is False

    def test_non_monotonic_x_is_not_deduplicated(self):
        """去重仅对比上一个 X 值，假设 X 单调递增。

        若 X 非单调（5→6→5），回退值会被当作新值接受，同一 X 值可能出现多个 Y。
        这不是设计保证——``try_accept_x_value`` 仅抑制**连续重复**的 X 值；
        "每个 X 值仅保留首个 Y" 的完整保证依赖调用方保证 custom X 单调递增。
        """
        r = DefinitionResolver()
        r.try_accept_x_value("acc", 5.0)
        r.try_accept_x_value("acc", 6.0)
        # 单调假设下不会出现 5→6→5；此处仅记录连续去重的实现行为
        assert r.try_accept_x_value("acc", 5.0) is True

    @pytest.mark.parametrize(
        "delta,expected",
        [(1e-9, False), (1e-7, True)],
        ids=["within-epsilon", "outside-epsilon"],
    )
    def test_epsilon_threshold(self, delta, expected):
        """epsilon（1e-8）内的微小差异视为相同，超出视为不同。"""
        r = DefinitionResolver()
        r.try_accept_x_value("loss", 1.0)
        assert r.try_accept_x_value("loss", 1.0 + delta) is expected

    def test_x_value_zero_accepted_then_rejected(self):
        """X=0.0 不受 falsy 影响。"""
        r = DefinitionResolver()
        assert r.try_accept_x_value("acc", 0.0) is True
        assert r.try_accept_x_value("acc", 0.0) is False

    def test_per_y_key_independent(self):
        """不同 Y key 各自独立去重。"""
        r = DefinitionResolver()
        r.try_accept_x_value("loss", 5.0)
        assert r.try_accept_x_value("acc", 5.0) is True  # acc 与 loss 独立


class TestStepSyncMergeReplace:
    """step_sync 在 define 时的默认值、merge、replace。

    默认值语义：从未显式指定时，step_sync = is_custom_x(effective.x_axis)，
    即提供 x_axis / step_metric 时默认 True，否则 False；显式 False 永远优先。
    """

    def test_default_step_sync_is_true(self):
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("loss", x_axis="epoch"))
        assert r.resolve_concrete("loss", "SCALAR").effective.step_sync is True

    def test_default_step_sync_false_without_x_axis(self):
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("loss", section_name="Train"))
        assert r.resolve_concrete("loss", "SCALAR").effective.step_sync is False

    def test_attaching_x_axis_later_defaults_step_sync_true(self):
        """先无 x_axis 定义、后补 x_axis 且未显式指定 step_sync → 默认 True。"""
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("loss", section_name="Train"))
        r.handle_define(MetricDefineEvent("loss", x_axis="epoch"))
        effective = r.resolve_concrete("loss", "SCALAR").effective
        assert effective.x_axis == "epoch"
        assert effective.step_sync is True

    def test_explicit_false_is_kept(self):
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("loss", x_axis="epoch", step_sync=False))
        assert r.resolve_concrete("loss", "SCALAR").effective.step_sync is False

    def test_merge_updates_step_sync_without_clearing_x_axis(self):
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("loss", x_axis="epoch"))
        r.handle_define(MetricDefineEvent("loss", step_sync=False))
        effective = r.resolve_concrete("loss", "SCALAR").effective
        assert effective.x_axis == "epoch"
        assert effective.step_sync is False

    def test_merge_omitted_step_sync_keeps_previous(self):
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("loss", x_axis="epoch", step_sync=False))
        r.handle_define(MetricDefineEvent("loss", section_name="Train"))
        effective = r.resolve_concrete("loss", "SCALAR").effective
        assert effective.step_sync is False
        assert effective.section_name == "Train"

    def test_overwrite_unspecified_resets_step_sync_to_true(self):
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("loss", x_axis="epoch", step_sync=False))
        r.handle_define(MetricDefineEvent("loss", x_axis="epoch", overwrite=True))
        assert r.resolve_concrete("loss", "SCALAR").effective.step_sync is True

    def test_overwrite_unspecified_resets_step_sync_to_false_without_x_axis(self):
        """overwrite 未指定 step_sync 时重置为默认值：无 custom X 则 False。"""
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("loss", x_axis="epoch"))
        r.handle_define(MetricDefineEvent("loss", overwrite=True))
        effective = r.resolve_concrete("loss", "SCALAR").effective
        assert effective.x_axis == "_step"
        assert effective.step_sync is False

    def test_overwrite_explicit_false_is_kept(self):
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("loss", x_axis="epoch"))
        r.handle_define(MetricDefineEvent("loss", x_axis="epoch", step_sync=False, overwrite=True))
        assert r.resolve_concrete("loss", "SCALAR").effective.step_sync is False

    def test_glob_step_sync_false_applies_to_matching_key(self):
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("train/*", is_glob=True, x_axis="epoch", step_sync=False))
        assert r.resolve_concrete("train/loss", "SCALAR").effective.step_sync is False

    def test_define_after_materialize_does_not_change_step_sync(self):
        r = DefinitionResolver()
        r.handle_define(MetricDefineEvent("loss", x_axis="epoch"))
        first = r.resolve_concrete("loss", "SCALAR")
        r.materialize_column("loss", "SCALAR", ColumnType.COLUMN_TYPE_SCALAR)
        r.handle_define(MetricDefineEvent("loss", step_sync=False))
        second = r.resolve_concrete("loss", "SCALAR")
        assert second is first
        assert second.effective.step_sync is True
