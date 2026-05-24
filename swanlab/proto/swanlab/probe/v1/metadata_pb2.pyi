from google.protobuf.internal import containers as _containers
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class MemoryUnit(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    MEMORY_UNIT_UNSPECIFIED: _ClassVar[MemoryUnit]
    MEMORY_UNIT_GB: _ClassVar[MemoryUnit]
    MEMORY_UNIT_MB: _ClassVar[MemoryUnit]

class AcceleratorVendor(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    ACCELERATOR_VENDOR_UNSPECIFIED: _ClassVar[AcceleratorVendor]
    ACCELERATOR_VENDOR_NVIDIA: _ClassVar[AcceleratorVendor]
    ACCELERATOR_VENDOR_ROCM: _ClassVar[AcceleratorVendor]
    ACCELERATOR_VENDOR_ILUVATAR: _ClassVar[AcceleratorVendor]
    ACCELERATOR_VENDOR_METAX: _ClassVar[AcceleratorVendor]
    ACCELERATOR_VENDOR_MOORETHREADS: _ClassVar[AcceleratorVendor]
    ACCELERATOR_VENDOR_ASCEND: _ClassVar[AcceleratorVendor]
    ACCELERATOR_VENDOR_CAMBRICON: _ClassVar[AcceleratorVendor]
    ACCELERATOR_VENDOR_HYGON: _ClassVar[AcceleratorVendor]
    ACCELERATOR_VENDOR_KUNLUNXIN: _ClassVar[AcceleratorVendor]
MEMORY_UNIT_UNSPECIFIED: MemoryUnit
MEMORY_UNIT_GB: MemoryUnit
MEMORY_UNIT_MB: MemoryUnit
ACCELERATOR_VENDOR_UNSPECIFIED: AcceleratorVendor
ACCELERATOR_VENDOR_NVIDIA: AcceleratorVendor
ACCELERATOR_VENDOR_ROCM: AcceleratorVendor
ACCELERATOR_VENDOR_ILUVATAR: AcceleratorVendor
ACCELERATOR_VENDOR_METAX: AcceleratorVendor
ACCELERATOR_VENDOR_MOORETHREADS: AcceleratorVendor
ACCELERATOR_VENDOR_ASCEND: AcceleratorVendor
ACCELERATOR_VENDOR_CAMBRICON: AcceleratorVendor
ACCELERATOR_VENDOR_HYGON: AcceleratorVendor
ACCELERATOR_VENDOR_KUNLUNXIN: AcceleratorVendor

class MetadataSnapshot(_message.Message):
    __slots__ = ("version", "hardware", "runtime", "git", "swanlab")
    VERSION_FIELD_NUMBER: _ClassVar[int]
    HARDWARE_FIELD_NUMBER: _ClassVar[int]
    RUNTIME_FIELD_NUMBER: _ClassVar[int]
    GIT_FIELD_NUMBER: _ClassVar[int]
    SWANLAB_FIELD_NUMBER: _ClassVar[int]
    version: int
    hardware: HardwareSnapshot
    runtime: RuntimeSnapshot
    git: GitSnapshot
    swanlab: SwanLabSnapshot
    def __init__(self, version: _Optional[int] = ..., hardware: _Optional[_Union[HardwareSnapshot, _Mapping]] = ..., runtime: _Optional[_Union[RuntimeSnapshot, _Mapping]] = ..., git: _Optional[_Union[GitSnapshot, _Mapping]] = ..., swanlab: _Optional[_Union[SwanLabSnapshot, _Mapping]] = ...) -> None: ...

class HardwareSnapshot(_message.Message):
    __slots__ = ("apple_silicon", "cpu", "memory", "accelerators")
    APPLE_SILICON_FIELD_NUMBER: _ClassVar[int]
    CPU_FIELD_NUMBER: _ClassVar[int]
    MEMORY_FIELD_NUMBER: _ClassVar[int]
    ACCELERATORS_FIELD_NUMBER: _ClassVar[int]
    apple_silicon: AppleSiliconSnapshot
    cpu: CPUSnapshot
    memory: MemorySnapshot
    accelerators: _containers.RepeatedCompositeFieldContainer[AcceleratorSnapshot]
    def __init__(self, apple_silicon: _Optional[_Union[AppleSiliconSnapshot, _Mapping]] = ..., cpu: _Optional[_Union[CPUSnapshot, _Mapping]] = ..., memory: _Optional[_Union[MemorySnapshot, _Mapping]] = ..., accelerators: _Optional[_Iterable[_Union[AcceleratorSnapshot, _Mapping]]] = ...) -> None: ...

class CPUSnapshot(_message.Message):
    __slots__ = ("brand", "physical_count", "logical_count")
    BRAND_FIELD_NUMBER: _ClassVar[int]
    PHYSICAL_COUNT_FIELD_NUMBER: _ClassVar[int]
    LOGICAL_COUNT_FIELD_NUMBER: _ClassVar[int]
    brand: str
    physical_count: int
    logical_count: int
    def __init__(self, brand: _Optional[str] = ..., physical_count: _Optional[int] = ..., logical_count: _Optional[int] = ...) -> None: ...

class MemorySnapshot(_message.Message):
    __slots__ = ("total", "total_unit")
    TOTAL_FIELD_NUMBER: _ClassVar[int]
    TOTAL_UNIT_FIELD_NUMBER: _ClassVar[int]
    total: int
    total_unit: MemoryUnit
    def __init__(self, total: _Optional[int] = ..., total_unit: _Optional[_Union[MemoryUnit, str]] = ...) -> None: ...

class AppleSiliconSnapshot(_message.Message):
    __slots__ = ("name", "memory", "memory_unit", "cpu_count")
    NAME_FIELD_NUMBER: _ClassVar[int]
    MEMORY_FIELD_NUMBER: _ClassVar[int]
    MEMORY_UNIT_FIELD_NUMBER: _ClassVar[int]
    CPU_COUNT_FIELD_NUMBER: _ClassVar[int]
    name: str
    memory: int
    memory_unit: MemoryUnit
    cpu_count: int
    def __init__(self, name: _Optional[str] = ..., memory: _Optional[int] = ..., memory_unit: _Optional[_Union[MemoryUnit, str]] = ..., cpu_count: _Optional[int] = ...) -> None: ...

class AcceleratorSnapshot(_message.Message):
    __slots__ = ("vendor", "devices", "version", "cuda_version", "cann_version")
    VENDOR_FIELD_NUMBER: _ClassVar[int]
    DEVICES_FIELD_NUMBER: _ClassVar[int]
    VERSION_FIELD_NUMBER: _ClassVar[int]
    CUDA_VERSION_FIELD_NUMBER: _ClassVar[int]
    CANN_VERSION_FIELD_NUMBER: _ClassVar[int]
    vendor: AcceleratorVendor
    devices: _containers.RepeatedCompositeFieldContainer[DeviceSnapshot]
    version: str
    cuda_version: str
    cann_version: str
    def __init__(self, vendor: _Optional[_Union[AcceleratorVendor, str]] = ..., devices: _Optional[_Iterable[_Union[DeviceSnapshot, _Mapping]]] = ..., version: _Optional[str] = ..., cuda_version: _Optional[str] = ..., cann_version: _Optional[str] = ...) -> None: ...

class DeviceSnapshot(_message.Message):
    __slots__ = ("index", "name", "memory", "memory_unit")
    INDEX_FIELD_NUMBER: _ClassVar[int]
    NAME_FIELD_NUMBER: _ClassVar[int]
    MEMORY_FIELD_NUMBER: _ClassVar[int]
    MEMORY_UNIT_FIELD_NUMBER: _ClassVar[int]
    index: int
    name: str
    memory: int
    memory_unit: MemoryUnit
    def __init__(self, index: _Optional[int] = ..., name: _Optional[str] = ..., memory: _Optional[int] = ..., memory_unit: _Optional[_Union[MemoryUnit, str]] = ...) -> None: ...

class RuntimeSnapshot(_message.Message):
    __slots__ = ("os", "os_pretty", "hostname", "pid", "cwd", "python_version", "python_verbose", "python_executable", "command")
    OS_FIELD_NUMBER: _ClassVar[int]
    OS_PRETTY_FIELD_NUMBER: _ClassVar[int]
    HOSTNAME_FIELD_NUMBER: _ClassVar[int]
    PID_FIELD_NUMBER: _ClassVar[int]
    CWD_FIELD_NUMBER: _ClassVar[int]
    PYTHON_VERSION_FIELD_NUMBER: _ClassVar[int]
    PYTHON_VERBOSE_FIELD_NUMBER: _ClassVar[int]
    PYTHON_EXECUTABLE_FIELD_NUMBER: _ClassVar[int]
    COMMAND_FIELD_NUMBER: _ClassVar[int]
    os: str
    os_pretty: str
    hostname: str
    pid: int
    cwd: str
    python_version: str
    python_verbose: str
    python_executable: str
    command: str
    def __init__(self, os: _Optional[str] = ..., os_pretty: _Optional[str] = ..., hostname: _Optional[str] = ..., pid: _Optional[int] = ..., cwd: _Optional[str] = ..., python_version: _Optional[str] = ..., python_verbose: _Optional[str] = ..., python_executable: _Optional[str] = ..., command: _Optional[str] = ...) -> None: ...

class GitSnapshot(_message.Message):
    __slots__ = ("remote_url", "branch", "commit")
    REMOTE_URL_FIELD_NUMBER: _ClassVar[int]
    BRANCH_FIELD_NUMBER: _ClassVar[int]
    COMMIT_FIELD_NUMBER: _ClassVar[int]
    remote_url: str
    branch: str
    commit: str
    def __init__(self, remote_url: _Optional[str] = ..., branch: _Optional[str] = ..., commit: _Optional[str] = ...) -> None: ...

class SwanLabSnapshot(_message.Message):
    __slots__ = ("version", "run_dir")
    VERSION_FIELD_NUMBER: _ClassVar[int]
    RUN_DIR_FIELD_NUMBER: _ClassVar[int]
    version: str
    run_dir: str
    def __init__(self, version: _Optional[str] = ..., run_dir: _Optional[str] = ...) -> None: ...
