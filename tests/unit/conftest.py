import os

import pytest

# 在 import swanlab 前清理所有 SWANLAB_ 环境变量，
# 避免模块级 Settings 常量受宿主环境影响。
# per-test 隔离由下方 autouse fixture 通过 monkeypatch 精确控制。
for _key in list(os.environ):
    if _key.startswith("SWANLAB_"):
        del os.environ[_key]


@pytest.fixture(autouse=True)
def isolate_sdk_environment(tmp_path, monkeypatch):
    """
    自动为每个测试用例隔离 SDK 全局环境，确保单个用例失败后不污染后续用例。

    Setup（测试前）：
      1. 隔离磁盘路径与环境变量
      2. 重置全局 Settings 覆盖

    Teardown（测试后，无论成功或失败均执行）：
      1. 清理 Run —— 需在 client 之前，online 模式的 finish() 依赖 client
      2. 清理 Client 单例
      3. 清理 logger
      4. 清理 callbacker
    """

    import swanlab
    from swanlab.sdk.internal.context.components.callbacker import global_callbacker
    from swanlab.sdk.internal.core_python import client
    from swanlab.sdk.internal.pkg import console
    from swanlab.sdk.internal.run import clear_run
    from swanlab.sdk.internal.run.components.config import reset as reset_config
    from swanlab.sdk.internal.settings import Settings, set_global_settings

    # 1. 路径与环境变量隔离
    # 清理可能干扰测试的常见环境变量
    for env_var in [
        "SWANLAB_API_KEY",
        "SWANLAB_API_HOST",
        "SWANLAB_WEB_HOST",
        "SWANLAB_ROOT",
        "SWANLAB_LOG_DIR",
    ]:
        monkeypatch.delenv(env_var, raising=False)

    # 切换当前工作目录，防止加载到项目根目录的 .env 或 swanlab.yaml
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path)

    # 2. 内存全局配置重置
    set_global_settings(Settings())

    # 3. 清理 Run 单例
    # 必须在 client 重置之前执行：online 模式的 run.finish() 可能需要 client 发送最后请求
    # 若 finish() 本身出错（如后台线程异常），直接强制清除引用，防止污染下一个用例
    if (run := swanlab.run) is not None:
        # noinspection PyBroadException
        try:
            run.finish()
        except Exception:  # noqa: E722
            clear_run()

    # 4. 清理 Client 单例
    if client.exists():
        client.reset()

    # 5. 清理 logger
    console.reset()

    # 6. 清理 config
    reset_config()

    # 7. 清理 callbacker
    for callback in global_callbacker.registered_callbacks:
        global_callbacker.remove_callback(callback.name)

    yield

    # 不将一些全局单例子放在teardown中，因为它们可能被monkeypatch覆盖
