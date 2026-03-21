

 上传线程本质上是个定时任务，需要一个缓冲区

 swanlab data 模块产出的是protobuf

 每次log产出若干protobuf，然后上传线程那块需要根据protobuf的类型聚合起来，统一发给house/server

规范要求: Type Hint 没有问题
---

## 已完成：uploader 模块重构

基于 protobuf Record 重构上传线程，实现在 `uploader/` 目录下，替代 `raw_uploader/` 的自定义 model 设计。

### 模块结构

- `model.py` — `UploadType` 枚举 + `classify_record()` 按 Record oneof 分类 + `FileModel` 聚合文件类 Record
- `batch.py` — `MetricDict` / `create_data()` / `trace_metrics()` 分批上传，通过 `get_client()` 获取客户端
- `upload.py` — 各类型上传函数（scalar/media/column/log/file）+ `UPLOAD_DISPATCH` 分发表
- `thread.py` — `ThreadPool` + `UploadCollector` + `RecordQueue` + `TimerFlag`，保留线程池防止数据丢失

### Record 类型映射

| UploadType | Record oneof |
|---|---|
| SCALAR_METRIC | metric (value=scalar) |
| MEDIA_METRIC | metric (value=images/audios/videos/texts/echarts) |
| COLUMN | column |
| LOG | console |
| FILE | metadata / requirements / conda / config |

### 快捷导入

`swanlab/core_python/__init__.py` 通过 `sys.modules` 别名注册，支持 `from swanlab.core_python.uploader import ...` 短路径。

### 修正：客户端调用与类型兼容

- `batch.py` — `trace_metrics()` 通过 `_get_client()` 获取全局单例 `Client`，响应使用 `resp.raw.status_code`（匹配 `ApiResponse` 契约）；`data` 参数类型修正为 `Optional[Union[MetricDict, dict, list]]`，兼容 Python 3.9
- `upload.py` — 移除未使用的 `client` 和 `MetricDict` 导入
- `test_upload.py` — mock 返回值从 `(None, MagicMock(status_code=200))` 改为 `MagicMock(raw=MagicMock(status_code=200))`，与 `ApiResponse` 一致
