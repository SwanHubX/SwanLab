import pytest

from swanlab.sdk.internal.context import callbacker
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import log
from swanlab.sdk.internal.run import clear_run, get_run, has_run
from swanlab.sdk.internal.settings import Settings, settings


@pytest.fixture(autouse=True)
def isolate_sdk_environment(tmp_path, monkeypatch):
    """
    自动为每个测试用例隔离 SDK 全局环境，确保单个用例失败后不污染后续用例。

    Setup（测试前）：
      1. 隔离磁盘路径与环境变量
      2. 重置 Settings 单例字段

    Teardown（测试后，无论成功或失败均执行）：
      1. 清理 SwanLabRun —— 需在 client 之前，cloud 模式的 finish() 依赖 client
      2. 清理 Client 单例
      3. 清理 logger
      4. 清理 callbacker
    """

    # ------------------------------------------------------------------ #
    # Setup
    # ------------------------------------------------------------------ #

    # 1. 路径与环境变量隔离
    # 清理可能干扰测试的常见环境变量
    for env_var in ["SWANLAB_API_KEY", "SWANLAB_API_HOST", "SWANLAB_WEB_HOST", "SWANLAB_ROOT", "SWANLAB_LOG_DIR"]:
        monkeypatch.delenv(env_var, raising=False)

    # 切换当前工作目录，防止加载到项目根目录的 .env 或 swanlab.yaml
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path)

    # 2. 内存单例重置 (Settings)
    # 创建一个全新的 Settings 实例，将其字段逐一写回全局单例
    # （不能直接替换指针，因为外部代码持有的是 `from ... import settings` 的直接引用）
    fresh_settings = Settings()
    for field_name in Settings.model_fields.keys():
        value = getattr(fresh_settings, field_name)
        object.__setattr__(settings, field_name, value)
    # 清除所有"字段已被修改"的 Pydantic 内部标记，让 settings 彻底回到初始状态
    object.__setattr__(settings, "__pydantic_fields_set__", set())
    object.__setattr__(settings, "__pydantic_extra__", None)

    # ------------------------------------------------------------------ #
    # 执行测试
    # ------------------------------------------------------------------ #
    yield

    # ------------------------------------------------------------------ #
    # Teardown（无论用例成功或失败均执行）
    # ------------------------------------------------------------------ #

    # 1. 清理 SwanLabRun 单例
    # 必须在 client 重置之前执行：cloud 模式的 run.finish() 可能需要 client 发送最后请求
    # 若 finish() 本身出错（如后台线程异常），直接强制清除引用，防止污染下一个用例
    if has_run():
        try:
            get_run().finish()
        except Exception:
            clear_run()

    # 2. 清理 Client 单例
    if client.exists():
        client.reset()

    # 3. 清理 logger
    log.reset()

    # 4. 清理 callbacker
    for callback in callbacker.registered_callbacks:
        callbacker.remove_callback(callback.name)
