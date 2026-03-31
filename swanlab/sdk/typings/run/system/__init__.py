"""
@author: cunyue
@file: __init__.py
@time: 2026/3/30 17:20
@description: SwanLab 系统信息相关类型定义
我们要求系统元数据是一个强类型约束的字典，否则在后续的查询和展示中会非常麻烦。
在命名方式上：
  - Snapshot → 描述硬件或系统状态的静态快照信息
"""

from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

# ──────────────────────────────────────────────
# 硬件快照子模型（静态信息，启动时采集一次）
# ──────────────────────────────────────────────


class CPUSnapshot(BaseModel):
    """CPU 静态信息"""

    brand: Optional[str] = None
    """CPU 品牌/型号，macOS 上为 None"""
    physical_count: Optional[int] = None
    """物理核心数"""
    logical_count: Optional[int] = None
    """逻辑核心数（含超线程）"""

    model_config = ConfigDict(frozen=True)


class MemorySnapshot(BaseModel):
    """内存静态信息"""

    total: Optional[int] = Field(default=None, gt=0)
    """总内存大小"""
    total_unit: Optional[Literal["GB", "MB"]] = None
    """内存单位"""

    @field_validator("total_unit", mode="before")
    @classmethod
    def normalize_unit(cls, v: Optional[str]) -> Optional[str]:
        return v.upper() if v else None

    @model_validator(mode="after")
    def validate_unit_dependency(self) -> "MemorySnapshot":
        if (self.total is None) != (self.total_unit is None):
            raise ValueError("total and total_unit must both be None or both be provided")
        return self

    model_config = ConfigDict(frozen=True)


class DeviceSnapshot(BaseModel):
    """单个加速器设备信息，GPU/NPU/MLU/DCU/XPU 均使用此模型"""

    index: Optional[int] = Field(default=None, ge=0)
    """设备序号"""
    name: Optional[str] = None
    """设备型号"""
    memory: Optional[int] = Field(default=None, ge=0)
    """显存/内存总量"""
    memory_unit: Optional[Literal["GB", "MB"]] = None
    """显存单位"""

    @field_validator("memory_unit", mode="before")
    @classmethod
    def normalize_unit(cls, v: Optional[str]) -> Optional[str]:
        return v.upper() if v else None

    @model_validator(mode="after")
    def validate_unit_dependency(self) -> "DeviceSnapshot":
        if (self.memory is None) != (self.memory_unit is None):
            raise ValueError("memory and memory_unit must both be None or both be provided")
        return self

    model_config = ConfigDict(frozen=True)


AcceleratorVendor = Literal[
    "nvidia",
    "rocm",
    "iluvatar",
    "metax",
    "moorethreads",
    "ascend",
    "cambricon",
    "hygon",
    "kunlunxin",
]
"""加速器厂商/驱动类型"""


class AcceleratorSnapshot(BaseModel):
    """通用加速器快照，GPU/NPU/MLU/DCU/XPU 共用此结构。

    vendor 特有字段用 Optional 扩展，不适用时保持 None。
    """

    vendor: AcceleratorVendor
    """加速器厂商/驱动类型"""
    devices: list[DeviceSnapshot] = Field(default_factory=list)
    """设备列表，每个物理设备一个条目"""
    version: Optional[str] = None
    """加速器驱动版本"""

    # vendor 特有字段
    cuda_version: Optional[str] = None
    """CUDA 版本号，仅 NVIDIA GPU 有效"""
    cann_version: Optional[str] = None
    """CANN 工具包版本，仅 Ascend NPU 有效"""

    model_config = ConfigDict(frozen=True)


class AppleSiliconSnapshot(BaseModel):
    """Apple Silicon 静态信息（CPU + 统一内存合并描述，不拆分）"""

    name: Optional[str] = None
    """芯片型号，如 'Apple M3 Pro'"""
    memory: Optional[int] = Field(default=None, gt=0)
    """统一内存总量"""
    memory_unit: Optional[Literal["GB", "MB"]] = None
    """内存单位"""
    cpu_count: Optional[int] = Field(default=None, gt=0)
    """CPU 核心数"""

    @field_validator("memory_unit", mode="before")
    @classmethod
    def normalize_unit(cls, v: Optional[str]) -> Optional[str]:
        return v.upper() if v else None

    @model_validator(mode="after")
    def validate_unit_dependency(self) -> "AppleSiliconSnapshot":
        if (self.memory is None) != (self.memory_unit is None):
            raise ValueError("memory and memory_unit must both be None or both be provided")
        return self

    model_config = ConfigDict(frozen=True)


class RuntimeSnapshot(BaseModel):
    """运行时环境静态信息"""

    os: Optional[str] = None
    """操作系统平台，如 'Linux-5.15.0'"""
    os_pretty: Optional[str] = None
    """发行版友好名称，如 'Ubuntu 22.04'"""
    hostname: Optional[str] = None
    pid: Optional[int] = None
    cwd: Optional[str] = None
    python_version: Optional[str] = None
    python_verbose: Optional[str] = None
    python_executable: Optional[str] = None
    command: Optional[str] = None
    """当前运行的命令行，包含所有参数"""

    model_config = ConfigDict(frozen=True)


class GitSnapshot(BaseModel):
    """Git 仓库信息"""

    remote_url: Optional[str] = None
    """远程仓库地址"""
    branch: Optional[str] = None
    commit: Optional[str] = None

    model_config = ConfigDict(frozen=True)


class HardwareSnapshot(BaseModel):
    """硬件相关的系统快照

    Apple Silicon 场景下 cpu/memory 为 None，apple_chip 有值。
    """

    apple_silicon: Optional[AppleSiliconSnapshot] = None
    cpu: Optional[CPUSnapshot] = None
    memory: Optional[MemorySnapshot] = None
    accelerators: list[AcceleratorSnapshot] = Field(default_factory=list)
    """所有加速器快照列表，支持异构场景"""

    model_config = ConfigDict(frozen=True)


# ──────────────────────────────────────────────
# MetadataSnapshot：完整的系统快照
# ──────────────────────────────────────────────


class MetadataSnapshot(BaseModel):
    """启动时采集的完整系统快照，用于上报展示。

    各字段为 None 表示未采集（被 Settings 关闭）或采集失败。
    """

    hardware: Optional[HardwareSnapshot] = None
    runtime: Optional[RuntimeSnapshot] = None
    git: Optional[GitSnapshot] = None

    def del_hardware(self):
        return self.model_copy(update={"hardware": None})

    model_config = ConfigDict(frozen=True)


# ──────────────────────────────────────────────
# SystemShim：monitor 工厂的唯一依赖
# ──────────────────────────────────────────────


class AcceleratorMonitorConfig(BaseModel):
    """单个加速器类型的监控配置"""

    vendor: AcceleratorVendor
    """加速器厂商/驱动类型，决定使用哪套采集接口"""
    device_indices: list[int]
    """需要监控的设备索引列表"""

    model_config = ConfigDict(frozen=True)


class SystemShim(BaseModel):
    """从 MetadataSnapshot 中提取的结构性事实。

    SystemShim 作为垫片，将 MetadataSnapshot 的完整信息转换为 monitor 工厂所需的简化配置

    这是 monitor/hardware/ 工厂函数的唯一输入，schema 固定，
    由 SwanLab 内部控制，用户无法修改。
    """

    slug: Literal["linux", "macos-intel", "macos-arm", "windows", "unknown"]
    """当前平台标识，目前主要和cpu memory监控实现相关，所以仅区分了macos平台下的apple chip和intel cpu"""
    enable_cpu: bool = False
    enable_memory: bool = False
    accelerators: list[AcceleratorMonitorConfig] = Field(default_factory=list)
    """各类加速器的监控配置"""

    model_config = ConfigDict(frozen=True)

    @classmethod
    def from_snapshot(cls, snapshot: MetadataSnapshot, platform: str) -> "SystemShim":
        """从 MetadataSnapshot 提取，不做额外的系统探测。"""
        # 1. 处理平台标识
        _slug: Literal["linux", "macos-intel", "macos-arm", "windows", "unknown"]
        if snapshot.hardware is None:
            _slug = "unknown"
        elif platform.startswith("linux"):
            _slug = "linux"
        elif platform.startswith("darwin"):
            _slug = "macos-arm" if snapshot.hardware.apple_silicon is not None else "macos-intel"
        else:
            _slug = "windows"

        # 2. 处理 CPU 和内存监控开关
        if snapshot.hardware is not None:
            enable_cpu = snapshot.hardware.cpu is not None or snapshot.hardware.apple_silicon is not None
            enable_memory = snapshot.hardware.memory is not None or snapshot.hardware.apple_silicon is not None
        else:
            enable_cpu = False
            enable_memory = False

        # 3. 处理加速器监控配置
        accelerators = []
        if snapshot.hardware is not None:
            for acc in snapshot.hardware.accelerators:
                if acc.vendor and acc.devices:
                    indices = [d.index for d in acc.devices if d.index is not None]
                    if indices:
                        accelerators.append(AcceleratorMonitorConfig(vendor=acc.vendor, device_indices=indices))

        return cls(
            enable_cpu=enable_cpu,
            enable_memory=enable_memory,
            accelerators=accelerators,
            slug=_slug,
        )


# ──────────────────────────────────────────────
# SystemInfo：完整容器
# ──────────────────────────────────────────────


class SystemEnvironment(BaseModel):
    """system/ 模块对外输出的完整信息容器。

    - shim:     monitor 工厂依赖，schema 固定
    - metadata: 系统元信息
    """

    shim: SystemShim
    metadata: MetadataSnapshot
    requirements: Optional[str]
    conda: Optional[str]
