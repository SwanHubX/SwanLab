from tutils import open_dev_mode
import swanlab

swanlab.login(open_dev_mode())

from fastai.vision.all import *
from swanlab.integration.fastai import SwanLabCallback


# 加载数据
path = untar_data(URLs.PETS)
dls = ImageDataLoaders.from_name_re(
    path, get_image_files(path / "images"), pat=r"([^/]+)_\d+.jpg$", item_tfms=Resize(224)
)

# 定义模型
learn = vision_learner(dls, resnet34, metrics=error_rate)

# 添加SwanLabCallback
learn.fit_one_cycle(
    5,
    cbs=SwanLabCallback(
        project="fastai-swanlab-integration-test",
        experiment_name="super-test",
        description="Test fastai integration with swanlab",
        logdir="./logs",
    ),
)
