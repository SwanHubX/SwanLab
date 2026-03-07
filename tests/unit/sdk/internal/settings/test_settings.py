"""
@author: cunyue
@file: test_settings.py
@time: 2026/3/5 15:56
@description: 测试 SwanLabSettings 基本方法
"""

from pathlib import Path

import pytest
from pydantic_settings import SettingsConfigDict

from swanlab.sdk.internal.settings import Settings


def test_path_validation(tmp_path):
    """测试默认值加载，以及路径是否会自动创建"""
    settings = Settings(root=Path(tmp_path) / "test", log_dir=Path(tmp_path) / "log")
    assert not settings.root.exists()
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


def test_mode_validation():
    """测试 mode 字段的校验与转换"""
    s_default = Settings()
    assert s_default.mode == "cloud"

    s_online = Settings(mode="online")  # type: ignore
    assert s_online.mode == "cloud"

    with pytest.raises(ValueError, match="Invalid mode: invalid, allowed values are"):
        Settings(mode="invalid")  # type: ignore


def test_url_resolution_logic(monkeypatch):
    """测试 api_host 和 web_host 的跨字段推导与路由/斜杠清理逻辑"""

    # 1. 默认情况：保持最新的默认值（无 /api 后缀）
    s_default = Settings()
    assert s_default.api_host == "https://api.swanlab.cn"
    assert s_default.web_host == "https://swanlab.cn"

    # 2. 仅自定义 api_host (带复杂路由)：自动清除路由，并顺便推导 web_host
    s_api = Settings(api_host="http://10.0.0.1:8080/api/v1/run/")
    assert s_api.api_host == "http://10.0.0.1:8080"
    assert s_api.web_host == "http://10.0.0.1:8080"

    # 3. 仅自定义 api_host (未带协议头)：自动补全 https://，清除路由，并推导 web_host
    s_api_no_scheme = Settings(api_host="custom-api.com/api/")
    assert s_api_no_scheme.api_host == "https://custom-api.com"
    assert s_api_no_scheme.web_host == "https://custom-api.com"

    # 4. 两者都自定义：清理各自格式，互不覆盖
    # api_host 会被清除路由，web_host 会被清理末尾斜杠
    s_both = Settings(web_host="http://web.local/", api_host="http://api.local/v1/")
    assert s_both.api_host == "http://api.local"
    assert s_both.web_host == "http://web.local"

    # 5. 仅自定义 web_host：清理斜杠，api_host 依然保持全局默认
    s_web = Settings(web_host="http://192.168.1.10/")
    assert s_web.web_host == "http://192.168.1.10"
    assert s_web.api_host == "https://api.swanlab.cn"


def test_url_env_resolution(monkeypatch):
    """测试环境变量注入时的 URL 路由清除与推导逻辑"""

    # 1. 仅通过环境变量注入 api_host (带复杂路由和斜杠)
    monkeypatch.delenv("SWANLAB_API_HOST", raising=False)
    monkeypatch.delenv("SWANLAB_WEB_HOST", raising=False)
    # 注意：如果你的 Settings 配置了 env_prefix="SWANLAB_"，这里应当是 SWANLAB_API_HOST
    monkeypatch.setenv("SWANLAB_API_HOST", "http://10.0.0.1:8080/api/v1/run/")

    s_env_api = Settings()
    # 验证：环境变量被正确读取，并且我们的 validate_urls 成功拦截并清洗了它
    assert s_env_api.api_host == "http://10.0.0.1:8080"
    assert s_env_api.web_host == "http://10.0.0.1:8080"

    # 清理环境变量，避免影响后续断言
    monkeypatch.delenv("SWANLAB_API_HOST", raising=False)

    # 2. 仅通过环境变量注入 web_host
    monkeypatch.setenv("SWANLAB_WEB_HOST", "http://env-web.local/")

    s_env_web = Settings()
    assert s_env_web.web_host == "http://env-web.local"
    assert s_env_web.api_host == "https://api.swanlab.cn"  # 保持默认不被推导

    monkeypatch.delenv("SWANLAB_WEB_HOST", raising=False)

    # 3. 测试优先级：代码显式传参 > 环境变量
    monkeypatch.setenv("SWANLAB_API_HOST", "http://env-api.local/api/")
    monkeypatch.setenv("SWANLAB_WEB_HOST", "http://env-web.local/")

    # 用户在代码里硬编码传入了 api_host
    s_mixed = Settings(api_host="http://code-api.local/api/v2/")

    # 验证：代码传参覆盖了 API_HOST 环境变量，且由于 API_HOST 的存在，
    # 我们配置的 WEB_HOST 环境变量被保留（清理了斜杠）
    assert s_mixed.api_host == "http://code-api.local"
    assert s_mixed.web_host == "http://env-web.local"
