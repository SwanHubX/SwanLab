"""
@author: caddiesnew
@file: __init__.py
@time: 2026/8/10
@description: define_metric resolver。

职责：
1. 管理 exact/glob rule 和 concrete state；
2. 在 MetricLogEvent 中解析 concrete definition、执行 step_sync 注入；
3. 为 definition 变化生成 ColumnRecord。
"""

from typing import Dict, Optional, Set, Tuple

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import (
    ColumnClass,
    ColumnRecord,
    ColumnType,
    SectionType,
)
from swanlab.sdk.internal.bus import MetricDefineEvent
from swanlab.sdk.internal.pkg import console, helper

from .state import (
    ConcreteState,
    DefinitionPatch,
    EffectiveDefinition,
    RuleState,
    make_effective,
)

__all__ = [
    "DefinitionResolver",
    "DefinitionPatch",
    "EffectiveDefinition",
    "RuleState",
    "ConcreteState",
]

# 每个 key 的 step 来源记录上限，防止长训练中 _x_step_origins 无限增长。
# 注入去重只关心"当前 step 是否已注入过"，历史 step 的来源记录无需长期保留。
_MAX_ORIGIN_STEPS_PER_KEY = 16

# custom X 值去重的 epsilon：|新值 - 上次值| < 此值视为相同。
_X_EPSILON = 1e-8


class DefinitionResolver:
    """define_metric resolver，生命周期与 Run 一致。

    由 BackgroundConsumer 单线程调用，无需加锁。
    """

    def __init__(self) -> None:
        # exact key → rule 状态（define_metric 中不含 '*' 的 key 注册于此）
        self._exact_rules: Dict[str, RuleState] = {}
        # glob pattern（末尾单 '*'）→ rule 状态，物化时按最长前缀匹配
        self._glob_rules: Dict[str, RuleState] = {}
        # (metric_class, key) → 首次 log 时物化的 concrete 定义快照；
        # first-writer-wins，之后的 define 不回溯修改已物化条目
        self._concrete: Dict[Tuple[str, str], ConcreteState] = {}
        # custom X 源 key → 最近一次真实 log 的 X 值，是跨 step 注入的取值来源
        self._custom_x_cache: Dict[str, float] = {}
        # custom X 源 key → {step → 来源标记 REAL / INJECTED / REAL_CONFLICT}；
        # 用于同 step 注入去重与"真实 X 晚于注入 X"的冲突判定（记录数有上限，见 _prune_origins）
        self._x_step_origins: Dict[str, Dict[int, str]] = {}
        # 已发过 X 冲突警告的 key，保证每 key 只警告一次
        self._warned_conflict_x: Set[str] = set()
        # 被某 rule 当作 custom X 轴的 key 集合。
        # 只有这些 key 才需要在 _custom_x_cache / _x_step_origins 中注册，
        # 避免对未参与 X 关系的海量 key 产生额外开销
        self._custom_x_keys: Set[str] = set()
        # Y key → 上次消费的 X 值；连续重复 X 值上的后续 Y 被去重（try_accept_x_value，epsilon 容差）
        self._y_last_x_value: Dict[str, float] = {}

    # ── rule 管理 ──────────────────────────────────────────────

    def handle_define(self, event: MetricDefineEvent) -> None:
        """处理 define_metric 事件，登记 exact/glob rule。

        glob 与 exact 都只登记 rule：图表定义遵循 first-writer-wins，由首次 log 时的
        materialize_column 决定；后续 define 不回溯已物化的 key、也不产出 ColumnRecord。
        同一 key 在 log 前多次 define 仍以最后一次为准（merge/replace 在 rule 上累积）。
        """
        # 1. 从事件提取 presence-aware 补丁，并按 is_glob 选定 rule 分区
        key = event.key
        patch = DefinitionPatch(
            x_axis=event.x_axis,
            section_name=event.section_name,
            hidden=event.hidden,
            step_sync=event.step_sync,
        )
        rules = self._glob_rules if event.is_glob else self._exact_rules

        # 2. 计算 effective：overwrite 从默认值 replace，否则在已有 rule（或默认）上 merge
        if event.overwrite:
            effective = self._replace_effective(patch)
        else:
            existing = rules.get(key)
            base = existing.effective if existing else self._default_effective()
            effective = self._merge_effective(base, patch)

        # 3. 登记 rule（exact / glob 分区，同 key 后一次 define 覆盖前一次）
        rules[key] = RuleState(
            identity=key,
            is_glob=event.is_glob,
            effective=effective,
            patch=patch,
            overwrite=event.overwrite,
        )

        # 4. 记录该 rule 引用的 custom X 源 key，供 update_custom_x_cache 按需登记
        if helper.is_custom_x(effective.x_axis):
            self._custom_x_keys.add(effective.x_axis)  # type: ignore[arg-type]

        # rule 登记完成。exact 不再立即重新物化已出现的 concrete：
        # 图表定义遵循 first-writer-wins，由首次 log 时的 materialize_column 决定。
        # 后续 define 只更新 rule，影响之后首次出现的 key。
        return

    @staticmethod
    def _default_effective() -> EffectiveDefinition:
        return make_effective("_step", None, False, True)

    @staticmethod
    def _merge_effective(base: EffectiveDefinition, patch: DefinitionPatch) -> EffectiveDefinition:
        x_axis = patch.x_axis if patch.x_axis is not None else base.x_axis
        section_name = patch.section_name if patch.section_name is not None else base.section_name
        hidden = patch.hidden if patch.hidden is not None else base.hidden
        step_sync = patch.step_sync if patch.step_sync is not None else base.step_sync
        return make_effective(x_axis, section_name, hidden, step_sync)

    @staticmethod
    def _replace_effective(patch: DefinitionPatch) -> EffectiveDefinition:
        x_axis = patch.x_axis if patch.x_axis is not None else "_step"
        section_name = patch.section_name
        hidden = patch.hidden if patch.hidden is not None else False
        step_sync = patch.step_sync if patch.step_sync is not None else True
        return make_effective(x_axis, section_name, hidden, step_sync)

    # ── concrete 解析 ──────────────────────────────────────────

    def resolve_concrete(self, key: str, metric_class: str) -> ConcreteState:
        """为 log 中的 key 解析 concrete definition。

        首次调用时按 exact → glob → automatic 优先级解析并注册 key；
        之后同一 (metric_class, key) 的所有调用直接返回首次快照，不再升级或回溯。
        define_metric 在 key 首次 log 后无法改变其定义
        """
        cache_key = (metric_class, key)
        existing = self._concrete.get(cache_key)
        if existing is not None:
            return existing

        # 以下仅在首次物化时执行

        # 1. exact 优先
        exact_rule = self._exact_rules.get(key)
        if exact_rule:
            effective = self._adjust_for_class(exact_rule.effective, metric_class, key)
            state = ConcreteState(
                key=key,
                metric_class=metric_class,
                effective=effective,
                source="exact",
                emitted_effective=None,
            )
            self._concrete[cache_key] = state
            return state

        # 2. 最长 prefix glob
        glob_rule = self._match_glob(key)
        if glob_rule:
            effective = self._adjust_for_class(glob_rule.effective, metric_class, key)
            state = ConcreteState(
                key=key,
                metric_class=metric_class,
                effective=effective,
                source="glob",
                emitted_effective=None,
            )
            self._concrete[cache_key] = state
            return state

        # 3. automatic 兜底
        effective = self._adjust_for_class(self._default_effective(), metric_class, key)
        state = ConcreteState(
            key=key,
            metric_class=metric_class,
            effective=effective,
            source="automatic",
            emitted_effective=None,
        )
        self._concrete[cache_key] = state
        return state

    def _match_glob(self, key: str) -> Optional[RuleState]:
        """找到最长 prefix glob rule。"""
        best: Optional[RuleState] = None
        best_len = -1
        for pattern, rule in self._glob_rules.items():
            prefix = pattern[:-1]
            if key.startswith(prefix) and len(prefix) > best_len:
                best = rule
                best_len = len(prefix)
        return best

    @staticmethod
    def _adjust_for_class(effective: EffectiveDefinition, metric_class: str, key: str) -> EffectiveDefinition:
        """根据 metric class 调整 effective definition。MEDIA 一律回退 _step；
        SCALAR 允许自引用（x_axis == key，图像为直线 y=x）。"""
        if metric_class == "MEDIA":
            return make_effective("_step", effective.section_name, effective.hidden, effective.step_sync)
        return effective

    # ── ColumnRecord 生成 ──────────────────────────────────────

    def materialize_column(
        self, key: str, metric_class: str, column_type: ColumnType, section_rule_index: int
    ) -> Optional[ColumnRecord]:
        """检查 effective 变化，返回 ColumnRecord 或 None。"""
        cache_key = (metric_class, key)
        state = self._concrete.get(cache_key)
        if state is None:
            return None
        state.column_type = int(column_type)
        if state.emitted_effective == state.effective:
            return None
        state.emitted_effective = state.effective
        return self._build_column_record(state, column_type, section_rule_index)

    @staticmethod
    def _build_column_record(state: ConcreteState, column_type: ColumnType, section_rule_index: int) -> ColumnRecord:
        effective = state.effective
        if state.source == "automatic":
            return ColumnRecord(
                column_class=ColumnClass.COLUMN_CLASS_CUSTOM,
                column_key=state.key,
                column_type=column_type,
                section_name=helper.calculate_section_name(state.key, section_rule_index),
                section_type=SectionType.SECTION_TYPE_PUBLIC,
            )
        col = ColumnRecord(
            column_class=ColumnClass.COLUMN_CLASS_CUSTOM,
            column_key=state.key,
            column_type=column_type,
            section_type=SectionType.SECTION_TYPE_PUBLIC,
        )
        col.section_name = effective.section_name or ""
        if state.metric_class == "SCALAR" and effective.x_axis != "_step":
            col.x_axis = effective.x_axis
        col.hidden = effective.hidden
        return col

    # ── custom X 值去重 ───────────────────────────────────────

    def try_accept_x_value(self, y_key: str, x_value: float) -> bool:
        """X 值与上次相同（epsilon 内）→ 重复返回 False；不同 → 记录并返回 True。

        类比 BaseMetric.try_accept_step(step)：同一 X 值上首次 Y 值为准，后续丢弃。

        .. note::
            仅抑制**连续重复**的 X 值。"每个 X 值仅保留首个 Y" 的完整保证依赖
            调用方保证 custom X 单调递增；非单调 X（如 5→6→5）下回退值会被当作新值
            接受，可能在同一 X 值上出现多个 Y 点。
        """
        last = self._y_last_x_value.get(y_key)
        if last is not None and abs(last - x_value) < _X_EPSILON:
            return False
        self._y_last_x_value[y_key] = x_value
        return True

    # ── step_sync: custom X cache ─────────────────────────────

    @property
    def has_custom_x(self) -> bool:
        """是否登记过任何 custom X 源 key。

        为 False 时消费端预扫描（transform + cache 更新）可整体跳过：
        update_custom_x_cache 本就以 _custom_x_keys 为门槛，explicit_scalars
        的消费点也全部位于 is_custom_x 分支之后。_custom_x_keys 仅由
        handle_define 添加、clear() 清空，FIFO 队列保证 define 先于后续 log。
        """
        return bool(self._custom_x_keys)

    def update_custom_x_cache(self, key: str, value: float, step: int) -> None:
        """缓存真实 custom X 值并登记 step 来源为 REAL。

        仅当 key 是某 rule 的 custom X 源（在 _custom_x_keys 中）时才登记，
        普通指标 key 无需占用 _custom_x_cache / _x_step_origins 内存。
        """
        if helper.is_system_key(key):
            return
        if key not in self._custom_x_keys:
            return
        self._custom_x_cache[key] = value
        origins = self._x_step_origins.setdefault(key, {})
        if step in origins and origins[step] == "INJECTED":
            # Y 先于 X 的场景：此 step 已有注入，真实 X 来了
            # 真实值仍进入 cache，但 Core 去重可能拒绝真实 X 落盘
            if key not in self._warned_conflict_x:
                console.warning(
                    f"Metric '{key}' at step {step}: real value arrives after injected value; "
                    f"the injected value at this step is kept. Log X before Y in the same step to avoid this."
                )
                self._warned_conflict_x.add(key)
            origins[step] = "REAL_CONFLICT"
        else:
            origins[step] = "REAL"
        self._prune_origins(origins)

    def get_custom_x(self, key: str) -> Optional[float]:
        return self._custom_x_cache.get(key)

    # ── step_sync: X 注入 ──────────────────────────────────────

    def try_inject_x(self, x_key: str, step: int) -> bool:
        """尝试为 (x_key, step) 登记注入。已有来源时返回 False。"""
        origins = self._x_step_origins.setdefault(x_key, {})
        if step in origins:
            return False
        origins[step] = "INJECTED"
        self._prune_origins(origins)
        return True

    @staticmethod
    def _prune_origins(origins: Dict[int, str]) -> None:
        """限制每个 key 保留的 step 来源记录数量，丢弃最旧的条目。

        注入去重只需判断"当前 step 是否已注入过"，历史 step 记录无长期价值；
        长训练中不做修剪会导致 _x_step_origins 随 step 数单调增长。
        """
        if len(origins) <= _MAX_ORIGIN_STEPS_PER_KEY:
            return
        # step 单调递增，按 key（step）升序丢弃最早的条目
        for stale_step in sorted(origins)[: len(origins) - _MAX_ORIGIN_STEPS_PER_KEY]:
            del origins[stale_step]

    # ── 清理 ───────────────────────────────────────────────────

    def clear(self) -> None:
        self._exact_rules.clear()
        self._glob_rules.clear()
        self._concrete.clear()
        self._custom_x_cache.clear()
        self._x_step_origins.clear()
        self._warned_conflict_x.clear()
        self._custom_x_keys.clear()
        self._y_last_x_value.clear()
