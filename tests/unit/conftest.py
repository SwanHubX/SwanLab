import pytest

from swanlab.sdk.internal.context import callbacker
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import log
from swanlab.sdk.internal.settings import Settings, settings


@pytest.fixture(autouse=True)
def isolate_sdk_environment(tmp_path, monkeypatch):
    """
    自动为每个测试用例隔离环境。
    1. 隔离磁盘路径 (SWANLAB_ROOT)
    2. 重置 Settings 单例
    3. 重置 Client 单例
    4. 清理 RunContext
    """

    # --- 1. 路径与环境变量隔离 ---
    # 清理可能干扰测试的常见环境变量
    for env_var in ["SWANLAB_API_KEY", "SWANLAB_API_HOST", "SWANLAB_WEB_HOST", "SWANLAB_ROOT", "SWANLAB_LOG_DIR"]:
        monkeypatch.delenv(env_var, raising=False)

    # 切换当前工作目录，防止加载到项目根目录的 .env 或 swanlab.yaml
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path)

    # --- 2. 内存单例重置 (Settings) ---
    # 创建一个全新的 Settings 实例
    fresh_settings = Settings()
    # 核心操作：将模块级的 settings 指针指向新实例
    # 测试结束后，monkeypatch 会自动将指针恢复原位
    # monkeypatch的做法只对settings.settings引用有效，但是我们一般是用"from swanlab.sdk.internal.settings import settings"
    # 这种引用是无效的，因此我们强制重制settings内部属性
    # monkeypatch.setattr(swanlab.sdk.internal.settings, "settings", fresh_settings)
    for field_name in Settings.model_fields.keys():
        # 直接从新实例中取出原始属性值（包括嵌套的子模型对象）
        # 这样塞回去的就是 ProjectSettings 对象，而不是 dict
        value = getattr(fresh_settings, field_name)
        object.__setattr__(settings, field_name, value)

        # 3. 强制重置 Pydantic 的内部追踪状态
        # 这一步能清除所有“字段已被修改”的标记，让 settings 彻底回到 Unset 状态
    object.__setattr__(settings, "__pydantic_fields_set__", set())
    # 清除私有属性缓存（如果有）
    object.__setattr__(settings, "__pydantic_extra__", None)

    # --- 3. 执行测试 ---
    yield

    # --- 4. 手动 Teardown 清理 ---
    # 确保测试结束后没有任何残留进入下一个用例
    # 1. 清理 Client 单例
    if client.exists():
        client.reset()
    # 2. 清理 logger
    log.reset()
    # 3. 清理 callbacker
    for callback in callbacker.registered_callbacks:
        callbacker.remove_callback(callback.name)
