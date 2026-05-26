from types import SimpleNamespace

import numpy as np

import swanlab.converter.mlf.sync as mlf_sync
import swanlab.converter.tfb.sync as tfb_sync
import swanlab.converter.wb.sync as wb_sync
from swanlab.converter.helper import extract_args


class _ConfigSpy:
    def __init__(self, updates):
        self.updates = updates

    def update(self, data=None, **kwargs):
        if data is not None:
            self.updates.append(data)
        if kwargs:
            self.updates.append(kwargs)


class _SwanLabSpy:
    def __init__(self):
        self.config_updates = []
        self.init_calls = []
        self.log_calls = []
        self.finished = 0
        self._has_run = False
        self.config = _ConfigSpy(self.config_updates)

    def has_run(self):
        return self._has_run

    def init(self, **kwargs):
        self._has_run = True
        self.init_calls.append(kwargs)

    def log(self, *args, **kwargs):
        self.log_calls.append((args, kwargs))

    def finish(self):
        self.finished += 1
        self._has_run = False


def test_extract_args_binds_positional_arguments_by_function_signature():
    def start_run(run_id=None, experiment_id=None, run_name=None):
        return run_id, experiment_id, run_name

    assert extract_args(
        start_run,
        ("run-123",),
        {"run_name": "actual-name"},
        ["run_name", "run_id", "experiment_id"],
    ) == ("actual-name", "run-123", None)


def test_sync_mlflow_accepts_keyword_experiment_id_and_extracts_run_name(monkeypatch):
    spy = _SwanLabSpy()

    class FakeMlflow:
        def set_experiment(self, experiment_name=None, experiment_id=None):
            return experiment_name, experiment_id

        def start_run(self, run_id=None, experiment_id=None, run_name=None):
            return run_id, experiment_id, run_name

        def end_run(self, status=None):
            return status

        def log_param(self, key, value):
            return key, value

        def log_params(self, params):
            return params

        def log_metric(self, key, value, step=None):
            return key, value, step

        def log_metrics(self, metrics, step=None):
            return metrics, step

    fake_mlflow = FakeMlflow()
    monkeypatch.setattr(mlf_sync, "vendor", SimpleNamespace(mlflow=fake_mlflow))
    monkeypatch.setattr(mlf_sync, "swanlab", spy)

    mlf_sync.sync_mlflow(mode="offline")

    fake_mlflow.set_experiment(experiment_id="exp-1")
    fake_mlflow.start_run("run-1", run_name="actual-name")

    assert spy.init_calls[-1]["name"] == "actual-name"


def test_sync_tensorboard_scalars_uses_real_parameter_names(monkeypatch):
    spy = _SwanLabSpy()

    class FakeSummaryWriter:
        def __init__(
            self,
            logdir=None,
            comment="",
            purge_step=None,
            max_queue=10,
            flush_secs=120,
            filename_suffix="",
            write_to_disk=True,
            log_dir=None,
            comet_config=None,
        ):
            self.logdir = logdir or log_dir

        def add_scalar(self, tag, scalar_value, global_step=None):
            return tag, scalar_value, global_step

        def add_scalars(self, main_tag, tag_scalar_dict, global_step=None):
            return main_tag, tag_scalar_dict, global_step

        def add_image(self, tag, img_tensor, global_step=None, walltime=None, dataformats="CHW"):
            return tag, img_tensor, global_step, walltime, dataformats

        def add_text(self, tag, text_string, global_step=None):
            return tag, text_string, global_step

        def close(self):
            return None

    monkeypatch.setattr(
        tfb_sync,
        "vendor",
        SimpleNamespace(tensorboardX=SimpleNamespace(SummaryWriter=FakeSummaryWriter)),
    )
    monkeypatch.setattr(tfb_sync, "swanlab", spy)

    tfb_sync.sync_tensorboardX()
    writer = FakeSummaryWriter(log_dir="runs/tb")
    writer.add_scalars(main_tag="train", tag_scalar_dict={"loss": 0.5}, global_step=7)

    assert spy.init_calls[-1]["config"] == {"tensorboard_logdir": "runs/tb"}
    assert spy.log_calls[-1][1] == {"data": {"train/loss": 0.5}, "step": 7}


def test_convert_tb_image_detaches_tensor_before_numpy(monkeypatch):
    class GradTensor:
        def __init__(self, detached=False):
            self.detached = detached

        def detach(self):
            return GradTensor(detached=True)

        def cpu(self):
            return self

        def numpy(self):
            if not self.detached:
                raise RuntimeError("requires grad")
            return np.ones((2, 3, 1))

    monkeypatch.setattr(tfb_sync, "vendor", SimpleNamespace(np=np))

    result = tfb_sync._convert_tb_image(GradTensor(), "HWC")

    assert result is not None
    assert result.shape == (2, 3, 1)


def test_sync_wandb_config_update_forwards_keyword_values(monkeypatch):
    spy = _SwanLabSpy()

    class FakeRun:
        def log(self, data=None, step=None, commit=None, sync=None):
            return data, step, commit, sync

        def finish(self):
            return None

    class FakeConfig:
        def update(self, d=None, allow_val_change=None, **kwargs):
            return d, allow_val_change, kwargs

    class FakeWandb:
        Image = object

        def __init__(self):
            self.sdk = SimpleNamespace(
                wandb_run=SimpleNamespace(Run=FakeRun),
                wandb_config=SimpleNamespace(Config=FakeConfig),
            )

        def init(
            self,
            entity=None,
            project=None,
            dir=None,
            id=None,
            name=None,
            notes=None,
            tags=None,
            config=None,
            config_exclude_keys=None,
            reinit=None,
            group=None,
            job_type=None,
        ):
            return entity, project, dir, id, name, notes, tags, config, config_exclude_keys, reinit, group, job_type

    fake_wandb = FakeWandb()
    monkeypatch.setattr(wb_sync, "vendor", SimpleNamespace(wandb=fake_wandb))
    monkeypatch.setattr(wb_sync, "swanlab", spy)

    wb_sync.sync_wandb()
    FakeConfig().update(lr=0.01, allow_val_change=True)

    assert {"lr": 0.01} in spy.config_updates
