from mmengine.config import Config
from mmengine.registry import VISUALIZERS
import mmengine
import swanlab

print(f"MMEngine Version: {mmengine.__version__}")
print(f"SwanLab Version: {swanlab.__version__}")

cfg_text = """
# swanlab visualizer
custom_imports = dict(  # 引入SwanLab作为日志记录器，对于部分不支持custom_imports的项目可以直接初始化SwanlabVisBackend并加入vis_backends
    imports=["swanlab.integration.mmengine"], allow_failed_imports=False
)

vis_backends = [
    dict(
        type="SwanlabVisBackend",
        init_kwargs={ # swanlab.init 参数
            "project": "swanlab-mmengine",
            "experiment_name": "Your exp",  # 实验名称
            "description": "Note whatever you want",  # 实验的描述信息
        },
    ),
]

visualizer = dict(
    type="Visualizer",
    vis_backends=vis_backends,
    name="visualizer",
)

"""

cfg = Config.fromstring(cfg_text, ".py")

custom_vis = VISUALIZERS.build(cfg.visualizer)
print(custom_vis)
