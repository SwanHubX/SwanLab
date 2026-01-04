"""
@author: Zhou QiYang
@file: run.py
@time: 2025/12/27 15:44
@description: 获取单个实验功能的单元测试
"""

import random
from unittest.mock import patch, MagicMock

import swanlab


# example_code:
#     api = swanlab.Api()
#     exp = api.run(path="username/project/expid")  # 可通过api.runs()获取expid
#     print(exp.__dict__)
#     print(exp.history())


def make_fake_single_experiment():
    return {
        "cuid": "".join(random.choices("0123456789abcdefghijklmnopqrstuvwxyz", k=20)),
        "name": "run-test-single",
        "createdAt": "2025-01-01T01:00:00Z",
        "description": "single experiment test",
    }


def make_fake_project_experiments():
    return [
        {
            "job": None,
            "cuid": "".join(random.choices("0123456789abcdefghijklmnopqrstuvwxyz", k=20)),
            "name": "run-test-single",
            "show": True,
            "pin": False,
            "state": "FINISHED",
            "colors": ["#b15fbb", "#b15fbb"],
            "cluster": "GG",
            "cloneType": "COPY",
            "rootProId": "".join(random.choices("0123456789abcdefghijklmnopqrstuvwxyz", k=10)),
            "rootExpId": "".join(random.choices("0123456789abcdefghijklmnopqrstuvwxyz", k=10)),
            "createdAt": "2025-01-01T01:00:00Z",
            "updatedAt": "2025-01-01T01:00:00Z",
            "finishedAt": "2025-01-01T03:00:00Z",
            "pinnedAt": None,
            "description": "single experiment test",
            "labels": [{"name": "label1"}, {"name": "label2"}],
            "user": {"username": "username", "name": "name"},
            "profile": {"config": {"lr": 0.00001, "loss": {"value": 123}}, "scalar": {"acc": 0.92, "loss": 0.5}},
            "resumes": [],
        }
    ]


def make_fake_history_csv(keys=None):
    if keys is None:
        keys = ['acc', 'loss']

    import pandas as pd

    data = dict()
    data['step'] = [i for i in range(100)]
    for key in keys:
        data[key] = [random.random() for _ in range(100)]
    data['timestamp'] = [f"2025-01-01T{i:02d}:00:00Z" for i in range(100)]

    df = pd.DataFrame(data)
    return df.to_csv(index=False).encode('utf-8')


def make_fake_metrics_url():
    return {"url": "https://fake-url.com/metrics.csv"}


def create_mock_response(csv_content):
    mock_response = MagicMock()
    mock_response.content = csv_content
    mock_response.__enter__ = MagicMock(return_value=mock_response)
    mock_response.__exit__ = MagicMock(return_value=None)
    return mock_response


class ApiExpContext:
    def __init__(self, fake_single_exp=None, fake_project_exps=None, fake_metrics_url=None):
        self.fake_single_exp = fake_single_exp or make_fake_single_experiment()
        self.fake_project_exps = make_fake_project_experiments() if fake_project_exps is None else fake_project_exps
        self.fake_metrics_url = fake_metrics_url
        self.patches = []
        self.api = None
        self.exp = None

    def __enter__(self):
        patch1 = patch("swanlab.api.api.get_single_experiment", return_value=self.fake_single_exp)
        patch2 = patch("swanlab.api.api.get_project_experiments", return_value=self.fake_project_exps)
        self.patches.extend([patch1, patch2])

        if self.fake_metrics_url is not None:
            patch3 = patch("swanlab.api.thread.get_experiment_metrics", return_value=self.fake_metrics_url)
            self.patches.append(patch3)

        for p in self.patches:
            p.start()

        self.api = swanlab.Api()
        self.exp = self.api.run(path="username/test-project/expid")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        for p in self.patches:
            p.stop()
        return False


def test_api_run_basic():
    """测试获取单个实验的基本属性完整性"""
    fake_project_exps = make_fake_project_experiments()
    fake_exp_data = fake_project_exps[0]

    with ApiExpContext(fake_project_exps=fake_project_exps) as ctx:
        exp = ctx.exp
        assert exp.name == fake_exp_data["name"]
        assert exp.state == fake_exp_data["state"]
        assert exp.description == fake_exp_data["description"]
        assert exp.group == fake_exp_data["cluster"]
        assert exp.job == fake_exp_data["job"]
        assert exp.created_at == fake_exp_data["createdAt"]

        assert exp.user.username == fake_exp_data["user"]["username"]
        assert exp.user.name == fake_exp_data["user"]["name"]

        fake_labels = fake_exp_data["labels"]
        assert len(exp.labels) == len(fake_labels)
        for i, label in enumerate(exp.labels):
            assert label.name == fake_labels[i]["name"]

        fake_profile = fake_exp_data["profile"]
        assert exp.config == fake_profile["config"]
        assert exp.summary == fake_profile["scalar"]

        fake_scalar = fake_profile["scalar"]
        expected_metric_keys = list(fake_scalar.keys())
        for key in expected_metric_keys:
            assert key in exp.metric_keys
        assert exp.history_line_count == 1
        assert exp.root_exp_id == fake_exp_data["rootExpId"]
        assert exp.root_pro_id == fake_exp_data["rootProId"]


def test_api_run_history():
    """测试实验history功能，包括获取所有指标、指定指标、x_axis、sample和pandas参数"""
    fake_metrics_url = make_fake_metrics_url()
    with ApiExpContext(fake_metrics_url=fake_metrics_url) as ctx:
        exp = ctx.exp
        fake_csv_acc = make_fake_history_csv(['acc'])
        fake_csv_loss = make_fake_history_csv(['loss'])
        call_count = [0]

        def mock_requests_get_all(url):
            call_count[0] += 1
            csv_content = fake_csv_acc if call_count[0] == 1 else fake_csv_loss
            return create_mock_response(csv_content)

        with patch("requests.get", side_effect=mock_requests_get_all):
            call_count[0] = 0
            history = exp.history()
            assert history is not None
            assert '_step' in history.columns
            assert 'acc' in history.columns
            assert 'loss' in history.columns
            assert len(history) == 100

        fake_csv_acc_only = make_fake_history_csv(['acc'])
        mock_response_acc = create_mock_response(fake_csv_acc_only)
        with patch("requests.get", return_value=mock_response_acc):
            history_keys = exp.history(keys=['acc'])
            assert history_keys is not None
            assert '_step' in history_keys.columns
            assert 'acc' in history_keys.columns
            assert 'loss' not in history_keys.columns

        call_count_x = [0]

        def mock_requests_get_x(url):
            call_count_x[0] += 1
            csv_content = fake_csv_acc if call_count_x[0] == 1 else fake_csv_loss
            return create_mock_response(csv_content)

        with patch("requests.get", side_effect=mock_requests_get_x):
            call_count_x[0] = 0
            history_x = exp.history(keys=['loss'], x_axis='acc')
            assert history_x is not None
            assert 'acc' in history_x.columns
            assert 'loss' in history_x.columns
            assert '_step' not in history_x.columns

        mock_response_sample = create_mock_response(fake_csv_acc_only)
        with patch("requests.get", return_value=mock_response_sample):
            history_sample = exp.history(keys=['acc'], sample=10)
            assert len(history_sample) == 10

        mock_response_dict = create_mock_response(fake_csv_acc_only)
        with patch("requests.get", return_value=mock_response_dict):
            history_dict = exp.history(keys=['acc'], pandas=False)
            assert isinstance(history_dict, list)
            assert len(history_dict) == 100
            if len(history_dict) > 0:
                assert isinstance(history_dict[0], dict)
                assert '_step' in history_dict[0] or 'step' in history_dict[0]
                assert 'acc' in history_dict[0]


def test_api_run_history_params():
    """测试history功能的sample、x_axis、pandas三个参数及其组合使用"""
    fake_metrics_url = make_fake_metrics_url()
    with ApiExpContext(fake_metrics_url=fake_metrics_url) as ctx:
        exp = ctx.exp
        fake_csv_acc = make_fake_history_csv(['acc'])
        fake_csv_loss = make_fake_history_csv(['loss'])

        mock_response_sample = create_mock_response(fake_csv_acc)
        with patch("requests.get", return_value=mock_response_sample):
            history_sample_10 = exp.history(keys=['acc'], sample=10)
            assert len(history_sample_10) == 10

            history_sample_50 = exp.history(keys=['acc'], sample=50)
            assert len(history_sample_50) == 50

            history_sample_none = exp.history(keys=['acc'], sample=None)
            assert len(history_sample_none) == 100

        call_count_x = [0]

        def mock_requests_get_x(url):
            call_count_x[0] += 1
            csv_content = fake_csv_acc if call_count_x[0] == 1 else fake_csv_loss
            return create_mock_response(csv_content)

        with patch("requests.get", side_effect=mock_requests_get_x):
            call_count_x[0] = 0
            history_x_axis = exp.history(keys=['loss'], x_axis='acc')
            assert history_x_axis is not None
            assert '_step' not in history_x_axis.columns
            assert history_x_axis.columns[0] == 'acc'
            assert 'loss' in history_x_axis.columns

        mock_response_pandas = create_mock_response(fake_csv_acc)
        with patch("requests.get", return_value=mock_response_pandas):
            import pandas as pd

            history_pandas_true = exp.history(keys=['acc'], pandas=True)
            assert isinstance(history_pandas_true, pd.DataFrame)
            assert '_step' in history_pandas_true.columns
            assert 'acc' in history_pandas_true.columns

            history_pandas_false = exp.history(keys=['acc'], pandas=False)
            assert isinstance(history_pandas_false, list)
            assert len(history_pandas_false) == 100
            if len(history_pandas_false) > 0:
                assert isinstance(history_pandas_false[0], dict)
                assert '_step' in history_pandas_false[0] or 'step' in history_pandas_false[0]
                assert 'acc' in history_pandas_false[0]

        call_count_combine = [0]

        def mock_requests_get_combine(url):
            call_count_combine[0] += 1
            csv_content = fake_csv_acc if call_count_combine[0] == 1 else fake_csv_loss
            return create_mock_response(csv_content)

        with patch("requests.get", side_effect=mock_requests_get_combine):
            call_count_combine[0] = 0
            history_combine = exp.history(keys=['loss'], x_axis='acc', sample=20, pandas=False)
            assert isinstance(history_combine, list)
            assert len(history_combine) == 20
            if len(history_combine) > 0:
                assert isinstance(history_combine[0], dict)
                assert 'acc' in history_combine[0]
                assert 'loss' in history_combine[0]
                assert '_step' not in history_combine[0]


def test_api_run_history_empty_keys():
    """测试没有指标时history返回空DataFrame"""
    fake_project_exps = make_fake_project_experiments()
    fake_exp = fake_project_exps[0]
    fake_exp["profile"]["scalar"] = {}

    with ApiExpContext(fake_project_exps=fake_project_exps) as ctx:
        exp = ctx.exp
        import pandas as pd

        history = exp.history()
        assert isinstance(history, pd.DataFrame)
        assert len(history) == 0
