// MetadataSnapshot 及其嵌套类型的 TypeScript 定义
//
// 对应 Python 定义：swanlab/sdk/typings/probe_python/__init__.py
// 序列化后的 JSON key 使用 Pydantic alias（如 _version）

// 加速器厂商/驱动类型
export type AcceleratorVendor =
  | "nvidia"       // NVIDIA GPU
  | "rocm"         // AMD ROCm
  | "iluvatar"     // 苍穹 Iluvatar
  | "metax"        // 壁仞 MetaX
  | "moorethreads" // 摩尔线程 Moore Threads
  | "ascend"       // 华为昇腾 Ascend NPU
  | "cambricon"    // 寒武纪 Cambricon MLU
  | "hygon"        // 海光 Hygon DCU
  | "kunlunxin";   // 昆仑芯 Kunlunxin XPU

// 单个加速器设备信息，GPU / NPU / MLU / DCU / XPU 均使用此模型
export interface DeviceSnapshot {
  index?: number | null;        // 设备序号
  name?: string | null;         // 设备型号
  memory?: number | null;       // 显存/内存总量
  memory_unit?: "GB" | "MB" | null; // 显存单位
}

// 通用加速器快照，GPU / NPU / MLU / DCU / XPU 共用此结构，vendor 特有字段用 Optional 扩展
export interface AcceleratorSnapshot {
  vendor: AcceleratorVendor;            // 加速器厂商/驱动类型
  devices: DeviceSnapshot[];            // 设备列表，每个物理设备一个条目
  version?: string | null;              // 加速器驱动版本
  cuda_version?: string | null;         // CUDA 版本号，仅 NVIDIA GPU 有效
  cann_version?: string | null;         // CANN 工具包版本，仅 Ascend NPU 有效
}

// CPU 静态信息
export interface CPUSnapshot {
  brand?: string | null;          // CPU 品牌/型号，macOS 上为 null
  physical_count?: number | null; // 物理核心数
  logical_count?: number | null;  // 逻辑核心数（含超线程）
}

// 内存静态信息
export interface MemorySnapshot {
  total?: number | null;                // 总内存大小
  total_unit?: "GB" | "MB" | null;      // 内存单位，与 total 必须同时存在或同时为 null
}

// Apple Silicon 静态信息，CPU + 统一内存合并描述，不拆分
// Apple Silicon 场景下 HardwareSnapshot.cpu/memory 为 null，apple_silicon 有值
export interface AppleSiliconSnapshot {
  name?: string | null;                 // 芯片型号，如 'Apple M3 Pro'
  memory?: number | null;               // 统一内存总量
  memory_unit?: "GB" | "MB" | null;     // 内存单位，与 memory 必须同时存在或同时为 null
  cpu_count?: number | null;            // CPU 核心数
}

// 硬件相关的系统快照
// Apple Silicon 场景下 cpu/memory 为 null，apple_silicon 有值
export interface HardwareSnapshot {
  apple_silicon?: AppleSiliconSnapshot | null; // Apple Silicon 芯片信息，非 Apple 设备为 null
  cpu?: CPUSnapshot | null;                    // CPU 信息
  memory?: MemorySnapshot | null;              // 内存信息
  accelerators: AcceleratorSnapshot[];         // 所有加速器快照列表，支持异构场景（如同时有 NVIDIA GPU + Ascend NPU）
}

// 运行时环境静态信息
export interface RuntimeSnapshot {
  os?: string | null;                // 操作系统平台，如 'Linux-5.15.0'
  os_pretty?: string | null;         // 发行版友好名称，如 'Ubuntu 22.04'
  hostname?: string | null;          // 主机名
  pid?: number | null;               // 进程 ID
  cwd?: string | null;               // 当前工作目录
  python_version?: string | null;    // Python 版本
  python_verbose?: string | null;    // Python 详细信息
  python_executable?: string | null; // Python 解释器路径
  command?: string | null;           // 当前运行的命令行，包含所有参数
}

// Git 仓库信息
export interface GitSnapshot {
  remote_url?: string | null; // 远程仓库地址
  branch?: string | null;     // 当前分支名
  commit?: string | null;     // 当前提交哈希
}

// 启动时采集的完整系统快照，用于上报展示
// 各字段为 null 表示未采集（被 Settings 关闭）或采集失败
export interface MetadataSnapshot {
  _version: number;                              // 快照 schema 版本号（序列化后 key 为 _version）
  hardware?: HardwareSnapshot | null;            // 硬件信息快照
  runtime?: RuntimeSnapshot | null;              // 运行时环境快照
  git?: GitSnapshot | null;                      // Git 仓库快照
}
