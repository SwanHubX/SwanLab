"""
@author: cunyue
@file: test_settings.py
@time: 2026/3/5 15:56
@description: 测试 SwanLabSettings 基本方法
"""

from pathlib import Path

from pydantic_settings import SettingsConfigDict

from swanlab.sdk.internal.settings import Settings


def test_path_validation():
    """测试默认值加载，以及路径是否会自动创建"""
    settings = Settings()
    assert settings.save_dir.exists()
    # log_dir 不会自动创建
    assert not settings.log_dir.exists()


def test_priority_yaml_over_env(tmp_path, monkeypatch):
    """测试优先级：本地 YAML > 环境变量"""

    # 设置低优先级的环境变量
    monkeypatch.setenv("SWANLAB_API_KEY", "env_key")
    monkeypatch.setenv("SWANLAB_HARDWARE_MONITOR", "false")

    s1 = Settings()
    assert s1.hardware.monitor is False
    assert s1.api_key == "env_key"

    # 创建高优先级的 swanlab.yaml
    yaml_content = "api_key: yaml_key\nlog_dir: ./custom_log\nhardware:\n  monitor: true"
    (tmp_path / "swanlab.yaml").write_text(yaml_content)

    s2 = Settings()
    assert s2.api_key == "yaml_key"
    assert s2.log_dir.name == "custom_log"
    assert s2.hardware.monitor is True


def test_merge_settings_dict_deep_update():
    """测试字典合并，确保深层嵌套不会被覆盖"""
    settings = Settings()

    # 仅修改 collect 下的 metadata，log_dir 应该保留
    settings.merge_settings({"hardware": {"monitor": False}})

    assert settings.hardware.monitor is False
    assert settings.hardware.interval == 10
    assert settings.log_dir.name == "swanlog"


def test_merge_settings_object_exclude_unset(tmp_path, monkeypatch):
    """测试对象合并，确保未设置的默认值不会覆盖已有配置"""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("SWANLAB_API_KEY", "env_key")
    settings = Settings(api_key="current_key")
    assert settings.api_key == "current_key"

    # 用户只传入了 log_dir，没有设置 metadata 和 api_key
    monkeypatch.delenv("SWANLAB_API_KEY")
    new_settings = Settings(log_dir=Path("./new_log"))
    settings.merge_settings(new_settings)

    assert settings.log_dir.name == "new_log"  # 新配置生效
    assert settings.hardware.monitor is True  # 原配置未被默认值 False 覆盖
    assert settings.api_key == "current_key"  # 原配置未被默认值 None 覆盖


def test_secrets_loading(tmp_path):
    """测试 K8s Secret 挂载文件读取"""
    secrets_dir = tmp_path / "secrets"
    secrets_dir.mkdir()
    (secrets_dir / "api_key").write_text("secret_from_k8s")

    # 通过继承临时覆写 secrets_dir 路径，避免修改全局环境变量
    class TestSettings(Settings):
        model_config = SettingsConfigDict(secrets_dir=str(secrets_dir))

    settings = TestSettings()
    assert settings.api_key == "secret_from_k8s"


def test_url_resolution_logic():
    """测试 api_url 和 web_url 的跨字段推导与斜杠清理逻辑"""

    # 1. 默认情况：都不改，保持默认
    s_default = Settings()
    assert s_default.web_url == "https://swanlab.cn"
    assert s_default.api_url == "https://api.swanlab.cn/api"

    # 2. 仅自定义 web_url：api_url 自动推导，并处理掉 web_url 末尾多余的斜杠
    s_web = Settings(web_url="http://10.0.0.1:8080/")
    assert s_web.web_url == "http://10.0.0.1:8080"
    assert s_web.api_url == "http://10.0.0.1:8080/api"

    # 3. 仅自定义 api_url：web_url 保持默认，清理 api_url 末尾斜杠
    s_api = Settings(api_url="http://custom-api.com/")
    assert s_api.web_url == "https://swanlab.cn"
    assert s_api.api_url == "http://custom-api.com"

    # 4. 两者都自定义：尊重用户的 api_url，不进行自动推导覆盖
    s_both = Settings(web_url="http://web.local", api_url="http://api.local/api/")
    assert s_both.web_url == "http://web.local"
    assert s_both.api_url == "http://api.local/api"
