"""
@author: cunyue
@file: test_init_e2e.py
@time: 2026/3/14
@description: 测试 swanlab.init() 的端到端行为

测试分组：
  - TestInitDisabledMode     : disabled 模式，无文件系统写入
  - TestInitLocalMode        : local 模式，验证目录创建
  - TestInitOfflineMode      : offline 模式，验证目录创建与 workspace 字段
  - TestInitReinit           : reinit 参数行为（跳过/重建）
  - TestInitSettingsPriority : 配置优先级（全局 < 自定义 < 传参）
  - TestInitResumeValidation : resume/id 校验逻辑
  - TestInitCloudMode        : cloud 模式，依赖本文件内的 HTTP mock fixtures
"""

import pytest
import responses as responses_lib

from swanlab.sdk.cmd.init import init
from swanlab.sdk.cmd.login import raw_login
from swanlab.sdk.cmd.merge_settings import merge_settings
from swanlab.sdk.internal.run import SwanLabRun, has_run
from swanlab.sdk.internal.settings import Settings

# ============================================================
# Cloud 模式测试常量
# ============================================================
API_HOST = "https://api.fake.swanlab.cn"
WEB_HOST = "https://test.swanlab.cn"
USERNAME = "test-user"
PROJECT = "test-project"
RUN_ID = "test-run-id"
API_KEY = "test-api-key"


# ============================================================
# Cloud 模式响应体数据工厂
# 返回与后端 schema 一致的 dict，支持 **overrides 局部定制
# ============================================================


def make_login_resp(**overrides) -> dict:
    """POST /api/login/api_key 响应体"""
    return {"sid": "mock-sid", "expiredAt": "2099-12-31T23:59:59.000Z", **overrides}


def make_init_project_resp(**overrides) -> dict:
    """POST /api/project 响应体（创建/获取项目）"""
    return {
        "name": PROJECT,
        "username": USERNAME,
        "path": f"/{USERNAME}/{PROJECT}",
        **overrides,
    }


def make_project_detail_resp(experiment_count: int = 0, **overrides) -> dict:
    """GET /api/project/{username}/{name} 响应体（项目详情）"""
    return {
        "name": PROJECT,
        "username": USERNAME,
        "path": f"/{USERNAME}/{PROJECT}",
        "visibility": "PRIVATE",
        "_count": {
            "experiments": experiment_count,
            "contributors": 1,
            "collaborators": 0,
            "clones": 0,
        },
        **overrides,
    }


def make_experiment_resp(**overrides) -> dict:
    """POST /api/project/{username}/{project}/experiment 响应体"""
    return {"cuid": RUN_ID, **overrides}


# ============================================================
# Cloud 模式 HTTP Mock Fixtures
#
# 分层设计：
#   层 1 - rsps：激活 responses.RequestsMock，所有 endpoint fixture 共享同一实例
#   层 2 - 单一 Endpoint Fixture：每个 fixture 只注册一个 API 端点
#   层 3 - 组合 Fixture：mock_cloud_init_apis 一次性装配 init() cloud 流程所需端点
#
# 未来 _init_cloud 新增 API 调用时，只需：
#   1. 新增 make_*_resp() 工厂函数
#   2. 新增对应 endpoint fixture
#   3. 将其追加为 mock_cloud_init_apis 的参数
# ============================================================


@pytest.fixture
def rsps():
    """激活 responses.RequestsMock 上下文，供所有 endpoint fixture 共享同一实例"""
    with responses_lib.RequestsMock() as r:
        yield r


@pytest.fixture
def mock_login_api(rsps):
    """注册 POST /api/login/api_key 端点"""
    rsps.add(responses_lib.POST, f"{API_HOST}/api/login/api_key", json=make_login_resp(), status=200)
    return rsps


@pytest.fixture
def mock_project_create_api(rsps):
    """注册 POST /api/project 端点（创建或获取项目）"""
    rsps.add(responses_lib.POST, f"{API_HOST}/api/project", json=make_init_project_resp(), status=201)
    return rsps


@pytest.fixture
def mock_project_get_api(rsps):
    """注册 GET /api/project/{username}/{name} 端点（项目详情）"""
    rsps.add(
        responses_lib.GET, f"{API_HOST}/api/project/{USERNAME}/{PROJECT}", json=make_project_detail_resp(), status=200
    )
    return rsps


@pytest.fixture
def mock_experiment_create_api(rsps):
    """注册 POST /api/project/{username}/{project}/experiment 端点（创建/恢复实验）"""
    rsps.add(
        responses_lib.POST,
        f"{API_HOST}/api/project/{USERNAME}/{PROJECT}/experiment",
        json=make_experiment_resp(),
        status=201,
    )
    return rsps


@pytest.fixture
def mock_cloud_settings():
    """将全局 settings 的 api_host/web_host 指向测试 HOST，避免意外触发生产环境"""
    merge_settings({"api_host": API_HOST, "web_host": WEB_HOST})


@pytest.fixture
def mock_cloud_init_apis(
    mock_cloud_settings,
    mock_login_api,
    mock_project_create_api,
    mock_project_get_api,
    mock_experiment_create_api,
    # 未来新增 API 端点 fixture 在此追加，例如：
    # mock_config_upload_api,
):
    """
    组合 fixture：一次性注册 init(mode='cloud') 当前所需的全部 HTTP 端点。
    所有子 fixture 通过 pytest 的 fixture 缓存机制共享同一个 rsps 实例。
    """
    pass


@pytest.fixture
def logged_in_client(mock_login_api, mock_cloud_settings):
    """
    调用 raw_login() 完成登录流程，确保全局 client 已创建。
    依赖 mock_login_api，故登录请求不会触及真实网络。
    清理由上级 conftest.py 的 isolate_sdk_environment 统一处理。
    """
    raw_login(api_key=API_KEY)
    yield


# ============================================================
# TestInitDisabledMode
# ============================================================


class TestInitDisabledMode:
    def test_init_returns_run_instance(self):
        """disabled 模式应返回 SwanLabRun 实例，且 has_run() 为 True"""
        run = init(mode="disabled")

        assert isinstance(run, SwanLabRun)
        assert has_run()

    def test_init_does_not_create_logdir(self):
        """disabled 模式不得创建日志目录"""
        run = init(mode="disabled")
        log_dir = run._ctx.config.settings.log_dir

        assert not log_dir.exists()

    def test_init_workspace_is_disabled(self):
        """disabled 模式下 workspace 字段应为 'disabled'"""
        run = init(mode="disabled")

        assert run._ctx.config.settings.project.workspace == "disabled"


# ============================================================
# TestInitLocalMode
# ============================================================


class TestInitLocalMode:
    def test_init_creates_logdir_with_gitignore(self):
        """local 模式应创建日志目录，并在首次创建时写入 .gitignore"""
        run = init(mode="local")
        log_dir = run._ctx.config.settings.log_dir

        assert log_dir.exists()
        assert (log_dir / ".gitignore").exists()

    def test_init_creates_run_subdirs(self):
        """local 模式应创建 run_dir、media_dir、files_dir、debug_dir"""
        run = init(mode="local")
        ctx = run._ctx

        assert ctx.run_dir.exists()
        assert ctx.media_dir.exists()
        assert ctx.files_dir.exists()
        assert ctx.debug_dir.exists()

    def test_init_workspace_is_local(self):
        """local 模式下 workspace 字段应为 'local'"""
        run = init(mode="local")

        assert run._ctx.config.settings.project.workspace == "local"

    def test_init_auto_generates_name_and_color(self):
        """local 模式未传 name/color 时，应自动生成非空字符串"""
        run = init(mode="local")
        s = run._ctx.config.settings

        assert s.experiment.name
        assert s.experiment.color


# ============================================================
# TestInitOfflineMode
# ============================================================


class TestInitOfflineMode:
    def test_init_creates_directories(self):
        """offline 模式应创建日志目录及运行子目录"""
        run = init(mode="offline")

        assert run._ctx.config.settings.log_dir.exists()
        assert run._ctx.run_dir.exists()

    def test_init_workspace_is_offline(self):
        """offline 模式下 workspace 字段应为 'offline'"""
        run = init(mode="offline")

        assert run._ctx.config.settings.project.workspace == "offline"


# ============================================================
# TestInitReinit
# ============================================================


class TestInitReinit:
    def test_reinit_false_returns_existing_run(self):
        """reinit=False（默认）时，第二次 init() 应返回已存在的 Run，而非新建"""
        run1 = init(mode="disabled")
        run2 = init(mode="disabled", reinit=False)

        assert run1 is run2

    def test_reinit_true_creates_new_run(self):
        """reinit=True 时，第二次 init() 应先 finish 旧 Run，再返回全新 Run"""
        run1 = init(mode="disabled")
        run2 = init(mode="disabled", reinit=True)

        assert run1 is not run2
        assert has_run()


# ============================================================
# TestInitSettingsPriority
# ============================================================


class TestInitSettingsPriority:
    def test_args_override_global_settings(self):
        """传参优先级最高：显式参数应覆盖通过 merge_settings 设置的全局配置"""
        merge_settings({"mode": "offline"})

        run = init(mode="disabled")

        assert run._ctx.config.settings.mode == "disabled"

    def test_custom_settings_override_global(self):
        """自定义 Settings 对象优先级高于全局配置，但低于显式传参"""
        merge_settings({"mode": "offline"})

        run = init(mode="disabled", settings=Settings(mode="local"))

        # 显式参数 mode="disabled" > custom mode="local" > global mode="offline"
        assert run._ctx.config.settings.mode == "disabled"


# ============================================================
# TestInitResumeValidation
# ============================================================


class TestInitResumeValidation:
    def test_resume_must_without_id_raises(self, monkeypatch):
        """cloud 模式 resume='must' 但未提供 id → AssertionError"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.prompt_init_mode", lambda _: ("cloud", True))

        with pytest.raises(AssertionError, match="Run id must be provided"):
            init(mode="cloud", resume="must")

    def test_resume_never_with_id_raises(self, monkeypatch):
        """cloud 模式 resume='never' 但同时提供了 id → AssertionError"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.prompt_init_mode", lambda _: ("cloud", True))

        with pytest.raises(AssertionError, match="Run id should not be provided"):
            init(mode="cloud", resume="never", id="some-run-id")

    def test_resume_validation_skipped_for_non_cloud_mode(self):
        """非 cloud 模式下，resume/id 校验不触发（即使传了也不报错）"""
        run = init(mode="offline", resume="must")

        assert isinstance(run, SwanLabRun)


# ============================================================
# TestInitCloudMode
# ============================================================


class TestInitCloudMode:
    def test_init_cloud_success(
        self, logged_in_client, mock_project_create_api, mock_project_get_api, mock_experiment_create_api
    ):
        """cloud 模式完整 init 流程：返回 SwanLabRun，has_run() 为 True"""
        run = init(mode="cloud", project=PROJECT)

        assert isinstance(run, SwanLabRun)
        assert has_run()

    def test_init_cloud_sets_project_and_workspace(
        self, logged_in_client, mock_project_create_api, mock_project_get_api, mock_experiment_create_api
    ):
        """cloud 模式下，workspace 和 project.name 应与后端响应同步"""
        run = init(mode="cloud", project=PROJECT)
        s = run._ctx.config.settings

        assert s.project.workspace == USERNAME
        assert s.project.name == PROJECT

    def test_init_cloud_project_already_exists(
        self, logged_in_client, rsps, mock_project_get_api, mock_experiment_create_api
    ):
        """POST /project 返回 409（项目已存在）时，应优雅降级为获取项目，继续初始化"""
        rsps.add(responses_lib.POST, f"{API_HOST}/api/project", json=make_init_project_resp(), status=409)

        run = init(mode="cloud", project=PROJECT)

        assert isinstance(run, SwanLabRun)
        assert run._ctx.config.settings.project.workspace == USERNAME

    def test_init_cloud_uses_all_expected_endpoints(
        self, logged_in_client, mock_project_create_api, mock_project_get_api, mock_experiment_create_api, rsps
    ):
        """验证 cloud init 确实调用了 project 和 experiment 端点"""
        init(mode="cloud", project=PROJECT)

        non_login_calls = [c.request.url for c in rsps.calls if "login" not in c.request.url]

        assert any("project" in u for u in non_login_calls), "应调用 project 端点"
        assert any("experiment" in u for u in non_login_calls), "应调用 experiment 端点"
