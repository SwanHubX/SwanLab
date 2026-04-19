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
  - TestInitFactoryDispatch  : 验证 factory 模式按模式分派组件类型
"""

import multiprocessing
from typing import cast

import pytest
import responses as responses_lib
import yaml

from swanlab.sdk.cmd.init import init
from swanlab.sdk.cmd.login import login_cli
from swanlab.sdk.cmd.merge_settings import merge_settings
from swanlab.sdk.internal.bus.emitter import RunEmitter
from swanlab.sdk.internal.run import Run, get_run, has_run
from swanlab.sdk.internal.run.config import config as global_config
from swanlab.sdk.internal.run.consumer import BackgroundConsumer
from swanlab.sdk.internal.run.factory import NullConsumer, NullEmitter
from swanlab.sdk.internal.settings import Settings, settings
from swanlab.sdk.typings.run import ModeType

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


def make_stop_experiment_resp(**overrides) -> dict:
    """PUT /api/project/{username}/{project}/runs/{cuid}/state 响应体"""
    return {"message": "ok", **overrides}


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
    """激活 responses.RequestsMock 上下文，供所有 endpoint fixture 共享同一实例。

    teardown 时在 responses 上下文关闭前自动 finish 当前 run，
    确保实验结束请求（stop_experiment）在 mock 活跃期间发出，
    避免 responses 报 "Not all requests have been executed" 错误。
    """
    with responses_lib.RequestsMock() as r:
        yield r
        # 在 responses 上下文关闭前，finish 当前 run 以触发 stop_experiment 请求
        if has_run():
            get_run().finish()


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
def mock_experiment_stop_api(rsps):
    """注册 PUT /api/project/{username}/{project}/runs/{cuid}/state 端点（停止实验）"""
    rsps.add(
        responses_lib.PUT,
        f"{API_HOST}/api/project/{USERNAME}/{PROJECT}/runs/{RUN_ID}/state",
        json=make_stop_experiment_resp(),
        status=200,
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
    mock_experiment_stop_api,
):
    """
    组合 fixture：一次性注册 init(mode='cloud') 当前所需的全部 HTTP 端点。
    所有子 fixture 通过 pytest 的 fixture 缓存机制共享同一个 rsps 实例。
    """
    pass


@pytest.fixture
def logged_in_client(mock_login_api, mock_cloud_settings):
    """
    调用 login_raw() 完成登录流程，确保全局 client 已创建。
    依赖 mock_login_api，故登录请求不会触及真实网络。
    清理由上级 conftest.py 的 isolate_sdk_environment 统一处理。
    """
    login_cli(api_key=API_KEY, host=settings.api_host)
    yield


# ============================================================
# TestInitDisabledMode
# ============================================================


class TestInitDisabledMode:
    def test_init_returns_run_instance(self):
        """disabled 模式应返回 Run 实例，且 has_run() 为 True"""
        run = init(mode="disabled")

        assert isinstance(run, Run)
        assert has_run()

    def test_init_does_not_create_logdir(self):
        """disabled 模式不得创建日志目录"""
        run = init(mode="disabled")
        log_dir = run._ctx.config.settings.log_dir

        assert not log_dir.exists()

    def test_init_workspace_is_none_by_default(self):
        """disabled 模式下 workspace 默认为 None"""
        run = init(mode="disabled")

        assert run._ctx.config.settings.project.workspace is None

    def test_init_custom_workspace_is_preserved(self):
        """disabled 模式下用户自定义 workspace 应被保留"""
        run = init(mode="disabled", workspace="my-workspace")

        assert run._ctx.config.settings.project.workspace == "my-workspace"

    def test_disabled_mode_does_not_write_media_files(self):
        """disabled 模式下 log 媒体数据不应写入任何本地文件"""
        import numpy as np

        run = init(mode="disabled")
        # 记录 log_dir 下所有文件（理论上不应存在任何文件）
        log_dir = run._ctx.config.settings.log_dir

        # log 各类媒体
        run.log({"loss": 0.5})
        run.log_image(key="img", data=np.zeros((10, 10, 3), dtype=np.uint8))
        run.log_text(key="txt", data="hello")
        run.log_audio(key="audio", data=np.zeros((1, 4410), dtype=np.float32), sample_rate=44100)
        run.finish()

        # log_dir 不应存在，或者即使存在也不应包含媒体文件
        assert not log_dir.exists()


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

    def test_init_workspace_is_none_by_default(self):
        """local 模式下 workspace 默认为 None"""
        run = init(mode="local")

        assert run._ctx.config.settings.project.workspace is None

    def test_init_custom_workspace_is_preserved(self):
        """local 模式下用户自定义 workspace 应被保留"""
        run = init(mode="local", workspace="my-team")

        assert run._ctx.config.settings.project.workspace == "my-team"

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

    def test_init_workspace_is_none_by_default(self):
        """offline 模式下 workspace 默认为 None"""
        run = init(mode="offline")

        assert run._ctx.config.settings.project.workspace is None

    def test_init_custom_workspace_is_preserved(self):
        """offline 模式下用户自定义 workspace 应被保留"""
        run = init(mode="offline", workspace="my-org")

        assert run._ctx.config.settings.project.workspace == "my-org"


# ============================================================
# TestInitReinit
# ============================================================


class TestInitReinit:
    def test_reinit_false_raises_when_run_exists(self):
        """reinit=False（默认）时，第二次 init() 应抛出 RuntimeError"""
        _ = init(mode="disabled")

        with pytest.raises(RuntimeError, match="`swanlab.init` requires an inactive Run"):
            init(mode="disabled", reinit=False)

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
        monkeypatch.setattr("swanlab.sdk.cmd.init.prompt_init_mode", lambda _: "cloud")

        with pytest.raises(AssertionError, match="Run id must be provided"):
            init(mode="cloud", resume="must")

    def test_resume_never_with_id_raises(self, monkeypatch):
        """cloud 模式 resume='never' 但同时提供了 id → AssertionError"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.prompt_init_mode", lambda _: "cloud")

        with pytest.raises(AssertionError, match="Run id should not be provided"):
            init(mode="cloud", resume="never", id="some-run-id")

    def test_resume_validation_skipped_for_non_cloud_mode(self):
        """非 cloud 模式下，resume/id 校验不触发（即使传了也不报错）"""
        run = init(mode="offline", resume="must")

        assert isinstance(run, Run)


# ============================================================
# TestInitCloudMode
# ============================================================


class TestInitCloudMode:
    def test_init_cloud_success(
        self,
        logged_in_client,
        mock_project_create_api,
        mock_project_get_api,
        mock_experiment_create_api,
        mock_experiment_stop_api,
    ):
        """cloud 模式完整 init 流程：返回 Run，has_run() 为 True"""
        run = init(mode="cloud", project=PROJECT)

        assert isinstance(run, Run)
        assert has_run()

    def test_init_cloud_sets_project_and_workspace(
        self,
        logged_in_client,
        mock_project_create_api,
        mock_project_get_api,
        mock_experiment_create_api,
        mock_experiment_stop_api,
    ):
        """cloud 模式下，workspace 和 project.name 应与后端响应同步"""
        run = init(mode="cloud", project=PROJECT)
        s = run._ctx.config.settings

        assert s.project.workspace == USERNAME
        assert s.project.name == PROJECT

    def test_init_cloud_project_already_exists(
        self, logged_in_client, rsps, mock_project_get_api, mock_experiment_create_api, mock_experiment_stop_api
    ):
        """POST /project 返回 409（项目已存在）时，应优雅降级为获取项目，继续初始化"""
        rsps.add(responses_lib.POST, f"{API_HOST}/api/project", json=make_init_project_resp(), status=409)

        run = init(mode="cloud", project=PROJECT)

        assert isinstance(run, Run)
        assert run._ctx.config.settings.project.workspace == USERNAME

    def test_init_cloud_uses_all_expected_endpoints(
        self,
        logged_in_client,
        mock_project_create_api,
        mock_project_get_api,
        mock_experiment_create_api,
        mock_experiment_stop_api,
        rsps,
    ):
        """验证 cloud init 确实调用了 project 和 experiment 端点"""
        init(mode="cloud", project=PROJECT)

        non_login_calls = [c.request.url for c in rsps.calls if "login" not in c.request.url]

        assert any("project" in u for u in non_login_calls), "应调用 project 端点"
        assert any("experiment" in u for u in non_login_calls), "应调用 experiment 端点"


# ============================================================
# TestInitConfigIntegration
# ============================================================


class TestInitConfigIntegration:
    def test_disabled_mode_does_not_bind_config(self):
        """disabled 模式下，config 不绑定，内存写入正常但不触发文件 IO"""
        global_config["lr"] = 0.01
        run = init(mode="disabled")

        assert run.config["lr"] == 0.01

    def test_local_mode_creates_config_file_after_init(self):
        """local 模式 init() 后，config.yaml 应已存在（即使 config 为空）"""
        run = init(mode="local")
        config_file = run._ctx.config_file

        assert config_file.exists()

    def test_pre_init_values_appear_in_config_file(self):
        """init() 前写入 config 的值，bindctx 时应一并 flush 到 config.yaml"""
        global_config["pre_lr"] = 0.001
        global_config["pre_epochs"] = 5

        run = init(mode="local")
        data = yaml.safe_load(run._ctx.config_file.read_text())

        assert data["pre_lr"]["value"] == 0.001
        assert data["pre_epochs"]["value"] == 5

    def test_post_init_values_written_to_config_file(self):
        """init() 后写入 config 的值，应实时写入 config.yaml"""
        run = init(mode="local")
        run.config["post_key"] = "hello"

        data = yaml.safe_load(run._ctx.config_file.read_text())
        assert data["post_key"]["value"] == "hello"

    def test_config_reset_after_finish(self):
        """finish() 后，全局 config 代理应指向 global config"""
        global_config["global_key"] = "global"
        run = init(mode="local")
        run.config["run_key"] = "run"
        run.finish()

        assert "global_key" in global_config
        assert "run_key" not in global_config

    def test_global_config_unchanged_after_finish(self):
        """finish() 后 global config 的键仍在"""
        global_config["persistent"] = "value"
        run = init(mode="local")
        run.finish()

        assert global_config["persistent"] == "value"

    def test_init_config_param_merged_into_run_config(self):
        """验证 init(config={"lr":0.01}) 写入 run.config"""
        run = init(mode="local", config={"lr": 0.01, "batch_size": 32})

        assert run.config["lr"] == 0.01
        assert run.config["batch_size"] == 32

    def test_proxy_points_to_run_config_after_init(self):
        """验证 global_config 在 run 活跃时指向 run.config"""
        run = init(mode="local")
        run.config["run_key"] = "run_value"

        assert global_config["run_key"] == "run_value"

    def test_offline_mode_creates_config_file(self):
        """offline 模式 init() 后，config.yaml 应已存在"""
        run = init(mode="offline")

        assert run._ctx.config_file.exists()

    def test_config_yaml_has_correct_structure(self):
        """config.yaml 的每个 key 应具有 {value, desc, sort} 结构"""
        global_config["lr"] = 0.01
        run = init(mode="local")

        data = yaml.safe_load(run._ctx.config_file.read_text())
        assert "value" in data["lr"]
        assert "desc" in data["lr"]
        assert "sort" in data["lr"]

    def test_reinit_rebinds_config(self):
        """reinit=True 重建 Run 后，config 应绑定到新的 config_file"""
        run1 = init(mode="local")
        run1.config["key1"] = "v1"
        config_file1 = run1._ctx.config_file

        run2 = init(mode="local", reinit=True)
        run2.config["key2"] = "v2"
        config_file2 = run2._ctx.config_file

        data2 = yaml.safe_load(config_file2.read_text())
        assert data2["key2"]["value"] == "v2"
        # 两次 run 的 config_file 路径不同（时间戳不同）
        assert config_file1 != config_file2


# ============================================================
# TestForkDetection
# ============================================================


@pytest.mark.skipif(not hasattr(__import__("os"), "register_at_fork"), reason="fork not available on this platform")
class TestForkDetection:
    """测试 fork 检测机制：register_at_fork 回调 + with_run 拦截 + has_run 失效"""

    def test_fork_sets_forked_flag(self):
        """真实 fork 后，子进程中 Run._forked 应为 True"""
        run = init(mode="disabled")
        assert not run._forked

        ctx = multiprocessing.get_context("fork")
        result = ctx.Queue()

        def child():
            result.put(run._forked)

        p = ctx.Process(target=child)
        p.start()
        p.join(timeout=5)
        assert result.get() is True

    def test_fork_clears_global_singleton(self):
        """真实 fork 后，子进程中 has_run() 应返回 False"""
        init(mode="disabled")
        assert has_run()

        ctx = multiprocessing.get_context("fork")
        result = ctx.Queue()

        def child():
            result.put(has_run())

        p = ctx.Process(target=child)
        p.start()
        p.join(timeout=5)
        assert result.get() is False

    def test_fork_log_raises_in_child(self):
        """真实 fork 后，子进程中 run.log() 应抛出 RuntimeError"""
        run = init(mode="disabled")

        ctx = multiprocessing.get_context("fork")
        result = ctx.Queue()

        def child():
            try:
                run.log({"loss": 0.5})
                result.put("no_error")
            except RuntimeError as e:
                result.put(str(e))

        p = ctx.Process(target=child)
        p.start()
        p.join(timeout=5)
        msg = result.get()
        assert "does not support fork" in msg

    def test_fork_finish_raises_in_child(self):
        """真实 fork 后，子进程中 finish() 也应抛出 RuntimeError"""
        run = init(mode="disabled")

        ctx = multiprocessing.get_context("fork")
        result = ctx.Queue()

        def child():
            try:
                run.finish()
                result.put("no_error")
            except RuntimeError as e:
                result.put(str(e))

        p = ctx.Process(target=child)
        p.start()
        p.join(timeout=5)
        msg = result.get()
        assert "does not support fork" in msg

    def test_fork_does_not_affect_parent(self):
        """真实 fork 后，父进程中的 Run 应不受影响"""
        run = init(mode="disabled")
        assert run.alive

        ctx = multiprocessing.get_context("fork")
        p = ctx.Process(target=lambda: None)
        p.start()
        p.join(timeout=5)

        # 父进程中 Run 仍然 alive
        assert run.alive
        assert not run._forked
        assert has_run()

    def test_can_init_new_run_in_child_after_fork(self):
        """真实 fork 后，子进程中可以重新 init 创建新 Run"""
        init(mode="disabled")

        ctx = multiprocessing.get_context("fork")
        result = ctx.Queue()

        def child():
            # fork 后 has_run() 为 False，可以重新 init
            try:
                run2 = init(mode="disabled")
                result.put(("ok", run2.alive, has_run()))
            except Exception as e:
                result.put(("error", str(e)))

        p = ctx.Process(target=child)
        p.start()
        p.join(timeout=5)
        r = result.get()
        assert r[0] == "ok"
        assert r[1] is True
        assert r[2] is True


# ============================================================
# TestEnsureRunDirE2E
# ============================================================


class TestEnsureRunDirE2E:
    """端到端测试：init 遇到 run_dir 冲突时能自动重试并成功初始化"""

    def test_init_local_retries_on_run_dir_conflict(self):
        """预先在 log_dir 下创建与当前时间戳同名的 run_dir，init 应自动重试并使用不同路径"""
        from datetime import datetime

        from swanlab.sdk.internal.settings import settings

        run_id = "conflict-test"
        log_dir = settings.log_dir
        log_dir.mkdir(parents=True, exist_ok=True)
        # 预先创建多个连续秒的冲突目录，确保命中
        all_conflict_names = []
        now = datetime.now()
        for offset in range(-1, 3):
            ts = datetime.fromtimestamp(now.timestamp() + offset)
            conflict_name = "run-" + ts.strftime("%Y%m%d_%H%M%S") + "-" + run_id
            (log_dir / conflict_name).mkdir(exist_ok=True)
            all_conflict_names.append(conflict_name)

        # init 应检测到冲突，等待后用新时间戳重试
        run = init(mode="offline", id=run_id)
        assert run._ctx.run_dir.exists()
        assert run._ctx.config.run_dir not in all_conflict_names
        run.finish()


# ============================================================
# TestInitFactoryDispatch
# 验证 factory 模式按模式分派组件类型的设计约定
# ============================================================


class TestInitFactoryDispatch:
    """验证 init() 创建的 Run 对象内部组件类型符合 factory 设计约定

    disabled 模式：NullEmitter + NullConsumer + unbound Config + 无 Monitor
    非 disabled 模式：RunEmitter + BackgroundConsumer + bound Config + 有 Monitor（环境允许时）
    """

    def test_disabled_uses_null_emitter(self):
        run = init(mode="disabled")
        assert isinstance(run._emitter, NullEmitter)

    @pytest.mark.parametrize("mode", ["local", "offline"])
    def test_non_disabled_uses_run_emitter(self, mode):
        run = init(mode=cast(ModeType, mode))
        assert isinstance(run._emitter, RunEmitter)

    def test_disabled_null_emitter_discards_events(self):
        """NullEmitter.emit() 丢弃事件，队列始终为空"""
        run = init(mode="disabled")
        run._emitter.emit("test-event")  # type: ignore
        assert run._emitter.queue.empty()

    def test_non_disabled_run_emitter_enqueues_events(self):
        """RunEmitter.emit() 将事件放入队列"""
        run = init(mode="local")
        run._emitter.emit("test-event")  # type: ignore
        assert not run._emitter.queue.empty()
        assert run._emitter.queue.get_nowait() == "test-event"

    def test_disabled_uses_null_consumer(self):
        run = init(mode="disabled")
        assert isinstance(run._consumer, NullConsumer)

    @pytest.mark.parametrize("mode", ["local", "offline"])
    def test_non_disabled_uses_background_consumer(self, mode):
        run = init(mode=cast(ModeType, mode))
        assert isinstance(run._consumer, BackgroundConsumer)

    def test_disabled_config_is_unbound(self):
        """disabled 模式 Config 不绑定文件，写操作仅内存"""
        run = init(mode="disabled")
        run.config["key"] = "value"

        assert run.config["key"] == "value"
        assert not run._ctx.config_file.exists()

    @pytest.mark.parametrize("mode", ["local", "offline"])
    def test_non_disabled_config_is_bound(self, mode):
        """非 disabled 模式 Config 绑定文件，写操作持久化"""
        run = init(mode=cast(ModeType, mode))
        run.config["key"] = "value"

        assert run._ctx.config_file.exists()

    def test_disabled_monitor_is_none(self):
        run = init(mode="disabled")
        assert run._monitor is None
