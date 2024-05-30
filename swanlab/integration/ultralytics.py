"""
Docs: https://docs.swanlab.cn/zh/guide_cloud/integration/integration-ultralytics.html

For adaptation to the ultralytics framework. Detailed usage are as follows:
------train.py in ultralytics------
from ultralytics import YOLO
from swanlab.integration.ultralytics import add_swanlab_callback

model = YOLO("yolov8n.pt")
add_swanlab_callback(model)

model.train(
    data="coco128.yaml",
    epoch=3,
)
---------------------------------

当使用Ultralytics做多卡、分布式训练时，需要在ultralytics的源码中添加回调代码。
在`your_env/ultralytics/utils/callbacks/base.py`的add_integration_callbacks函数中添加三行代码：

------.../ultralytics/utils/callbacks/base.py------
def add_integration_callbacks(instance):
    ...

    if "Trainer" in instance.__class__.__name__:
        ...
        from swanlab.integration.ultralytics import return_swanlab_callback

        sw_cb = return_swanlab_callback()

        callbacks_list.extend([ ..., sw_cb])
        ...
---------------------------------

"""

from ultralytics.models import YOLO
from ultralytics.utils.torch_utils import model_info_for_loggers
from collections import Counter
from typing import Optional, Dict, Any
import swanlab


_processed_plots = {}


class UltralyticsSwanlabCallback:
    def __init__(
        self,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        experiment_name: Optional[str] = None,
        description: Optional[str] = None,
        logdir: Optional[str] = None,
        mode: Optional[bool] = None,
        **kwargs: Any,
    ) -> None:
        self.step_counter = Counter()
        self._run = None

        self._swanlab_init: Dict[str, Any] = {
            "project": project,
            "workspace": workspace,
            "experiment_name": experiment_name,
            "description": description,
            "logdir": logdir,
            "mode": mode,
        }

        self._swanlab_init.update(**kwargs)

        self._project = self._swanlab_init.get("project")
        self._workspace = self._swanlab_init.get("workspace")
        self._experiment_name = self._swanlab_init.get("experiment_name")
        self._description = self._swanlab_init.get("decsription")
        self._logdir = self._swanlab_init.get("logdir")
        self._mode = self._swanlab_init.get("mode")

    def _log_plots(self, plots: dict, step: int, tag: str):
        """记录指标绘图和推理图像"""
        image_list = []
        for (
            name,
            params,
        ) in plots.copy().items():
            timestamp = params["timestamp"]
            if _processed_plots.get(name) != timestamp:
                image_list.append(swanlab.Image(str(name), caption=name.stem))
                _processed_plots[name] = timestamp

        if image_list:
            self._run.log({tag: image_list}, step=step)

    def on_pretrain_routine_start(self, trainer):
        """初始化实验记录器"""
        if swanlab.get_run() is None:
            self._run = swanlab.init(
                project=trainer.args.project if self._project is None else self._project,
                workspace=self._workspace,
                experiment_name=trainer.args.name if self._experiment_name is None else self._experiment_name,
                config=vars(trainer.args),
                description=self._description,
                logdir=self._logdir,
                mode=self._mode,
            )
        else:
            self._run = swanlab.get_run()
            self._run.config.update(vars(trainer.args))

    def on_train_epoch_end(self, trainer):
        """每个train epoch结束记录指标（仅训练）"""
        swanlab.log(trainer.label_loss_items(trainer.tloss, prefix="train"), step=trainer.epoch + 1)
        swanlab.log(trainer.lr, step=trainer.epoch + 1)

        if trainer.epoch == 1:
            self._log_plots(trainer.plots, step=trainer.epoch + 1, tag="Plots")

    def on_fit_epoch_end(self, trainer):
        """每个epoch结束记录指标和绘图（含训练和验证）"""
        self.step_counter[trainer.epoch + 1] += 1

        # final_eval阶段
        if self.step_counter[trainer.epoch + 1] == 2:
            for key, value in trainer.metrics.items():
                swanlab.log({f"FinalValBestModel/{key}": value}, step=trainer.epoch + 1)
        # train阶段
        else:
            swanlab.log(trainer.metrics, step=trainer.epoch + 1)
            if trainer.epoch + 1 == trainer.epochs:
                self.final_eval = True

        self._log_plots(trainer.plots, step=trainer.epoch + 1, tag="Train/Plots")
        self._log_plots(trainer.validator.plots, step=trainer.epoch + 1, tag="Train/ValPlots")

        if trainer.epoch == 0:
            self._run.log(model_info_for_loggers(trainer), step=trainer.epoch + 1)

    def on_train_end(self, trainer):
        """结束训练"""
        self._log_plots(trainer.plots, step=trainer.epoch + 1, tag="TrainEnd/Plots")
        self._log_plots(trainer.validator.plots, step=trainer.epoch + 1, tag="TrainEnd/ValPlots")
        self._run.finish()


def add_swanlab_callback(
    model: YOLO,
    project: Optional[str] = None,
    workspace: Optional[str] = None,
    experiment_name: Optional[str] = None,
    description: Optional[str] = None,
    logdir: Optional[str] = None,
    mode: Optional[bool] = None,
):
    ultralytics_swanlabcallback = UltralyticsSwanlabCallback(
        project=project,
        workspace=workspace,
        experiment_name=experiment_name,
        description=description,
        logdir=logdir,
        mode=mode,
    )

    """给Ultralytics模型添加swanlab回调函数"""
    callbacks = {
        "on_pretrain_routine_start": ultralytics_swanlabcallback.on_pretrain_routine_start,
        "on_fit_epoch_end": ultralytics_swanlabcallback.on_fit_epoch_end,
        "on_train_epoch_end": ultralytics_swanlabcallback.on_train_epoch_end,
        "on_train_end": ultralytics_swanlabcallback.on_train_end,
    }

    for event, callback_fn in callbacks.items():
        model.add_callback(event, callback_fn)

    return model


def return_swanlab_callback(
    project: Optional[str] = None,
    workspace: Optional[str] = None,
    experiment_name: Optional[str] = None,
    description: Optional[str] = None,
    logdir: Optional[str] = None,
    mode: Optional[bool] = None,
):
    ultralytics_swanlabcallback = UltralyticsSwanlabCallback(
        project=project,
        workspace=workspace,
        experiment_name=experiment_name,
        description=description,
        logdir=logdir,
        mode=mode,
    )

    callbacks = (
        {
            "on_pretrain_routine_start": ultralytics_swanlabcallback.on_pretrain_routine_start,
            "on_train_epoch_end": ultralytics_swanlabcallback.on_train_epoch_end,
            "on_fit_epoch_end": ultralytics_swanlabcallback.on_fit_epoch_end,
            "on_train_end": ultralytics_swanlabcallback.on_train_end,
        }
        if swanlab
        else {}
    )

    return callbacks
