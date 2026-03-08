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

用户可以通过merge_settings动态合并配置，但是在设计上，在执行`swanlab.init`和`swanlab.finish`之间，无法使用merge_settings。
"""

import os
from pathlib import Path
from typing import Any, ClassVar, Dict, Optional, Tuple, Type, Union, get_args
from urllib.parse import urlparse

from pydantic import Field, field_validator
from pydantic.functional_validators import model_validator
from pydantic_settings import (
    BaseSettings,
    PydanticBaseSettingsSource,
    SecretsSettingsSource,
    SettingsConfigDict,
    YamlConfigSettingsSource,
)

from swanlab.sdk.typings.run import ModeType

from .experiment import ExperimentSettings, ProjectSettings, RunSettings
from .integration import IntegrationSettings
from .metadata import ConsoleSettings, EnvSettings, HardwareSettings

__all__ = ["Settings", "settings", "strip_none"]

# 根据环境变量自动设置 secrets_dir
# 如果强制设置，会出现警告：https://github.com/pydantic/pydantic/issues/2175
secrets_dir_env = os.getenv("SWANLAB_SECRETS_DIR")
SECRETS_DIR: Optional[str] = secrets_dir_env or None


def strip_none(data: dict) -> dict:
    """
    递归剔除字典中的 None 值和空字典，主要用于配置项合并时的空值处理
    """
    clean_data = {}
    for k, v in data.items():
        if isinstance(v, dict):
            cleaned_v = strip_none(v)
            # 只有当嵌套字典里真的有非 None 的值时，才保留这个 key
            if cleaned_v:
                clean_data[k] = cleaned_v
        elif v is not None:
            clean_data[k] = v
    return clean_data


def root_factory() -> Path:
    # 向下兼容旧版本环境变量
    return Path(os.environ.get("SWANLAB_SAVE_DIR", str(Path.home() / ".swanlab")))


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

    interactive: bool = True
    """
    Whether to enable interactive mode.
    If False, all user input prompts and related interactions will be disabled.
    Useful for CI/CD environments or background batch jobs.
    """

    mode: ModeType = "cloud"
    """
    SwanLab Run mode.

    * `local`: Run SwanLab locally.
    * `cloud`: Run SwanLab on the cloud.
    * `disabled`: Disable SwanLab.
    * `offline`: Run SwanLab in offline mode.
    """

    @field_validator("mode", mode="before")
    def validate_mode(cls, v: Any) -> ModeType:
        if v in list(get_args(ModeType)):
            return v
        if v == "online":
            return "cloud"
        raise ValueError(f"Invalid mode: {v}, allowed values are {list(get_args(ModeType))}")

    root: Path = Field(default_factory=root_factory)
    """
    Directory for SwanLab saved files.
    """

    # @field_validator("save_dir", mode="before")
    # def validate_save_dir(cls, v: Union[str, Path]) -> Path:
    #     """在 Pydantic 校验它是 DirectoryPath 之前，先把它建出来"""
    #     path_v = Path(v)
    #
    #     if not path_v.exists():
    #         path_v.mkdir(parents=True, exist_ok=True)
    #     return path_v

    log_dir: Path = Field(default=Path.cwd() / "swanlog", validate_default=True)
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
    api_host: str = Field(default="https://api.swanlab.cn")
    """
    Base URL for SwanLab API services.
    """
    web_host: str = Field(default="https://swanlab.cn")
    """
    Base URL for SwanLab web services.
    It just a display URL for SwanLab web services, no actual effect on SDK behavior.
    """

    @model_validator(mode="before")
    @classmethod
    def strip_non_empty(cls, data: Dict) -> Dict:
        """
        删除空值和空字典，以适配传入None的情况，一般情况下此校验必须在其他model_validator之前定义
        如果出现部分字段需要识别None值，则在此校验之前定义model_validator
        """
        if isinstance(data, dict):
            data = strip_none(data)
        return data

    @model_validator(mode="before")
    @classmethod
    def validate_hosts(cls, data: Dict) -> Dict:
        """
        校验并清理 HOST 字段，确保它们以正确的格式存在
        在设计上，api_host 是最基础URL，但是有时候展示的前端URL和后端URL可能不一致
        所以在处理时，我们优先使用 api_host，然后根据需要（当没有显式配置 web_host 时）推导 web_host
        """
        if isinstance(data, dict):
            if "api_host" in data and data["api_host"]:
                raw_api = str(data["api_host"]).strip().rstrip("/")

                # 1. 如果没有携带协议头，默认拼接 https://
                if not raw_api.startswith(("http://", "https://")):
                    raw_api = f"https://{raw_api}"

                # 2. 交给 urlparse 解析，此时必然有 scheme 和 netloc
                parsed = urlparse(raw_api)
                # 重新拼接，完美清除掉所有的 path/query 等冗余信息
                clean_api = f"{parsed.scheme}://{parsed.netloc}"

                # 将纯净的 URL 写回
                data["api_host"] = clean_api

                # 3. 当且仅当没有显式配置 web_host 时，自动推导 web_host
                if "web_host" not in data:
                    data["web_host"] = clean_api

            # 统一处理 web_host 末尾的斜杠
            current_web: str = data.get("web_host", cls.model_fields["web_host"].default)
            if current_web and current_web.endswith("/"):
                data["web_host"] = current_web.rstrip("/")

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
