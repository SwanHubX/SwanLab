"""
@author: cunyue
@file: __init__.py
@time: 2026/3/5 14:38
@description: SwanLab 包配置项，根据优先级从低到高加载配置：
1. 默认值
2. 环境变量
3. 当前目录下 .env 文件
4. /etc/swanlab/*.{yaml,yml}
5. 当前目录下 swanlab.{yaml,yml}
6. K8S/Docker 容器 Secret 配置项文件

在设计上 Settings 仅是与用户交互的配置入口，不包含业务逻辑，这意味着仅检查必要的类型和格式和必要的默认值，不产生副作用：
1. 文件夹创建
2. 具体业务逻辑，如实验id生成与格式校验、实验名称长度校验等
"""

import os
from pathlib import Path
from typing import ClassVar, Dict, Optional, Tuple, Type, Union

from pydantic import DirectoryPath, Field, field_validator
from pydantic.functional_validators import model_validator
from pydantic_settings import (
    BaseSettings,
    PydanticBaseSettingsSource,
    SecretsSettingsSource,
    SettingsConfigDict,
    YamlConfigSettingsSource,
)

from swanlab.sdk.types.run import ModeType

from .experiment import ExperimentSettings, ProjectSettings, RunSettings
from .integration import IntegrationSettings
from .metadata import ConsoleSettings, EnvSettings, HardwareSettings

__all__ = ["Settings", "settings"]

# 根据环境变量自动设置 secrets_dir
# 如果强制设置，会出现警告：https://github.com/pydantic/pydantic/issues/2175
secrets_dir_env = os.getenv("SWANLAB_SECRETS_DIR")
SECRETS_DIR: Optional[str] = secrets_dir_env or None


class Settings(BaseSettings):
    Project: ClassVar[Type[ProjectSettings]] = ProjectSettings
    Run: ClassVar[Type[RunSettings]] = RunSettings
    Experiment: ClassVar[Type[ExperimentSettings]] = ExperimentSettings
    Hardware: ClassVar[Type[HardwareSettings]] = HardwareSettings
    Console: ClassVar[Type[ConsoleSettings]] = ConsoleSettings
    Env: ClassVar[Type[EnvSettings]] = EnvSettings
    Integration: ClassVar[Type[IntegrationSettings]] = IntegrationSettings

    debug: bool = False
    """
    Whether to enable debug mode for SwanLab.
    If enabled, SwanLab will output more detailed logs.
    """

    mode: ModeType = "cloud"
    """
    SwanLab Run mode.

    * `local`: Run SwanLab locally.
    * `cloud`: Run SwanLab on the cloud.
    * `disabled`: Disable SwanLab.
    * `offline`: Run SwanLab in offline mode.
    """

    save_dir: DirectoryPath = Field(default=Path.home() / ".swanlab", validate_default=True)
    """
    Directory for SwanLab saved files.
    """

    @field_validator("save_dir", mode="before")
    def validate_save_dir(cls, v: Union[str, Path]) -> Path:
        """在 Pydantic 校验它是 DirectoryPath 之前，先把它建出来"""
        path_v = Path(v)

        if not path_v.exists():
            path_v.mkdir(parents=True, exist_ok=True)
        return path_v

    log_dir: Path = Field(default=Path("./swanlog"), validate_default=True)
    """
    Directory for SwanLab logs.
    Semantically, this is just a path representation and the directory may NOT exist when loaded. 
    The actual folder creation is deferred to the SDK initialization phase to avoid side effects.
    """

    # 不在此处创建文件夹，而是在 sdk init 时创建
    # @field_validator("log_dir", mode="before")
    # def validate_log_dir(cls, v: Union[str, Path]) -> Path:
    #     """在 Pydantic 校验它是 DirectoryPath 之前，先把它建出来"""
    #     path_v = Path(v)
    #     # 即使它是一个已存在的同名文件，这里不报错，交给后面的 DirectoryPath 去报错
    #     if not path_v.exists():
    #         path_v.mkdir(parents=True, exist_ok=True)
    #     return path_v

    api_key: Optional[str] = Field(default=None)
    """
    API key for SwanLab services.
    """
    api_url: str = Field(default="https://api.swanlab.cn/api")
    """
    Base URL for SwanLab API services.
    """
    web_url: str = Field(default="https://swanlab.cn")
    """
    Base URL for SwanLab web services.
    If you set `web_url`, `api_url` will be automatically derived from it.
    """

    @model_validator(mode="before")
    @classmethod
    def validate_urls(cls, data: Dict) -> Dict:
        if isinstance(data, dict):
            web_url: str = data.get("web_url", cls.model_fields["web_url"].default)
            # 1. 当且仅当用户显式配置了 web_url，而没有显式配置 api_url
            # 注意：在 before 模式下，我们要检查 data 这个字典里是否存在对应的 key
            if "web_url" in data:
                clean_web = web_url.rstrip("/")
                data["web_url"] = clean_web
                if "api_url" not in data:
                    data["api_url"] = f"{clean_web}/api"

            # 2. 统一处理 api_url 末尾的斜杠
            # 即使是默认值，我们也需要处理它（如果 data 中没有，则取默认值处理）
            current_api: str = data.get("api_url", cls.model_fields["api_url"].default)
            if current_api and current_api.endswith("/"):
                data["api_url"] = current_api.rstrip("/")

        return data

    project: ProjectSettings = Field(default_factory=ProjectSettings)
    """
    Configuration for the project of this SwanLab run.
    """
    experiment: ExperimentSettings = Field(default_factory=ExperimentSettings)
    """
    Configuration for the experiment of this SwanLab run.
    """
    run: RunSettings = Field(default_factory=RunSettings)
    """
    Configuration for the run of this SwanLab experiment.
    """
    hardware: HardwareSettings = Field(default_factory=HardwareSettings)
    """
    Configuration for SwanLab hardware monitoring.
    """
    env: EnvSettings = Field(default_factory=EnvSettings)
    """
    Configuration for SwanLab environment information collection.
    """
    console: ConsoleSettings = Field(default_factory=ConsoleSettings)
    """
    Configuration for SwanLab terminal log collection.
    """
    integration: IntegrationSettings = Field(default_factory=IntegrationSettings)
    """
    Configuration for SwanLab integrations, including webhook, dashboard, etc.
    """

    model_config = SettingsConfigDict(
        env_prefix="SWANLAB_",
        env_nested_delimiter="_",
        env_file=".env",
        env_file_encoding="utf-8",
        env_nested_max_split=1,
        extra="ignore",
        frozen=True,
        # 指定 Secret 文件存放目录，通常在容器中挂载到这里
        secrets_dir=SECRETS_DIR,
        validate_assignment=True,
        str_strip_whitespace=True,
    )

    @classmethod
    def settings_customise_sources(
        cls,
        settings_cls: Type[BaseSettings],
        init_settings: PydanticBaseSettingsSource,
        env_settings: PydanticBaseSettingsSource,
        dotenv_settings: PydanticBaseSettingsSource,
        file_secret_settings: PydanticBaseSettingsSource,
    ) -> Tuple[PydanticBaseSettingsSource, ...]:

        # 优先级由高到低排列体现在返回的sources顺序：
        # 1. init_settings (merge_settings 传入的参数)
        # 2. 当前目录下 swanlab.yaml
        # 3. /etc/swanlab/*.yaml
        # 4. .env 文件
        # 5. file_secret_settings (容器 Secrets)
        # 6. env_settings (环境变量)
        # 7. 默认值 (Model Default)
        sources = [init_settings]

        # 5. 当前目录下 swanlab.{yaml,yml}
        for ext in ["yaml", "yml"]:
            local_file = Path(f"swanlab.{ext}")
            if local_file.exists():
                sources.append(YamlConfigSettingsSource(settings_cls, yaml_file=local_file))

        # 4. /etc/swanlab/*.{yaml,yml}
        etc_dir = Path("/etc/swanlab")
        if etc_dir.exists() and etc_dir.is_dir():
            etc_files = sorted(list(etc_dir.glob("*.yaml")) + list(etc_dir.glob("*.yml")), reverse=True)
            for file in etc_files:
                sources.append(YamlConfigSettingsSource(settings_cls, yaml_file=file))

        # 3. .env 文件
        sources.append(dotenv_settings)

        # file_secret_settings (Secrets 文件)
        # 对于Secrets文件，不需要额外的前缀，直接使用默认的环境变量前缀
        secrets_dir = settings_cls.model_config.get("secrets_dir")
        if secrets_dir:
            custom_secret_source = SecretsSettingsSource(
                settings_cls,
                secrets_dir=secrets_dir,
                env_prefix="",
            )
            sources.append(custom_secret_source)
        # 优先级高于普通环境变量，防止敏感信息被低优先级的 Env 覆盖
        sources.append(file_secret_settings)

        # 2. 环境变量
        sources.append(env_settings)

        return tuple(sources)

    def merge_settings(self, other: Union["Settings", dict]) -> None:
        """
        合并自定义设置
        """
        if isinstance(other, self.__class__):
            # 1. 使用 exclude_unset=True 提取用户显式设置的字段
            # 这样可以确保 SwanLabSettings("log_dir"="...") 不会带上默认的值
            update_data = other.model_dump(exclude_unset=True)
        elif isinstance(other, dict):
            update_data = other
        else:
            raise TypeError(f"Only {self.__class__.__name__} or dict can be merged, not {type(other)}")

        # 2. 获取当前实例的完整状态字典
        current_data = self.model_dump()

        # 3. 递归合并数据 (确保嵌套的 dict 如 collect.metadata 不会被整个替换)
        merged_data = _deep_update(current_data, update_data)

        # 4. 验证新数据：通过类构造函数生成临时实例以触发校验和 Path 转换
        # 这一步能确保传入的路径字符串被 field_validator 处理成 Path 对象并创建目录
        validated_instance = self.__class__(**merged_data)

        # 5. 更新当前单例
        for field_name in self.__class__.model_fields.keys():
            new_value = getattr(validated_instance, field_name)
            # 绕过 frozen=True 的限制
            object.__setattr__(self, field_name, new_value)

        # 由于我们 绕过了 frozen=True 的限制，需要手动同步 __pydantic_fields_set__
        object.__setattr__(self, "__pydantic_fields_set__", validated_instance.__pydantic_fields_set__)


def _deep_update(base_dict: dict, update_dict: dict) -> dict:
    """递归合并字典，用于嵌套模型的 merge_settings"""
    for k, v in update_dict.items():
        if k in base_dict and isinstance(base_dict[k], dict) and isinstance(v, dict):
            base_dict[k] = _deep_update(base_dict[k], v)
        else:
            base_dict[k] = v
    return base_dict


settings = Settings()
