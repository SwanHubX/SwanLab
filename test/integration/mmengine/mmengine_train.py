import torch.nn.functional as F
import torchvision
import torchvision.transforms as transforms
from torch.optim import SGD
from torch.utils.data import DataLoader

from mmengine.evaluator import BaseMetric
from mmengine.model import BaseModel
from mmengine.runner import Runner
from mmengine.visualization import Visualizer


from tutils import open_dev_mode
import swanlab
from swanlab.integration.mmengine import SwanlabVisBackend

swanlab.login(open_dev_mode())


class MMResNet50(BaseModel):
    def __init__(self):
        super().__init__()
        self.resnet = torchvision.models.resnet50()

    def forward(self, imgs, labels, mode):
        x = self.resnet(imgs)
        if mode == "loss":
            return {"loss": F.cross_entropy(x, labels)}
        elif mode == "predict":
            return x, labels


class Accuracy(BaseMetric):
    def process(self, data_batch, data_samples):
        score, gt = data_samples
        self.results.append(
            {
                "batch_size": len(gt),
                "correct": (score.argmax(dim=1) == gt).sum().cpu(),
            }
        )

    def compute_metrics(self, results):
        total_correct = sum(item["correct"] for item in results)
        total_size = sum(item["batch_size"] for item in results)
        return dict(accuracy=100 * total_correct / total_size)


norm_cfg = dict(mean=[0.491, 0.482, 0.447], std=[0.202, 0.199, 0.201])
train_dataloader = DataLoader(
    batch_size=32,
    shuffle=True,
    dataset=torchvision.datasets.CIFAR10(
        "data/cifar10",
        train=True,
        download=True,
        transform=transforms.Compose(
            [
                transforms.RandomCrop(32, padding=4),
                transforms.RandomHorizontalFlip(),
                transforms.ToTensor(),
                transforms.Normalize(**norm_cfg),
            ]
        ),
    ),
)

val_dataloader = DataLoader(
    batch_size=32,
    shuffle=False,
    dataset=torchvision.datasets.CIFAR10(
        "data/cifar10",
        train=False,
        download=True,
        transform=transforms.Compose([transforms.ToTensor(), transforms.Normalize(**norm_cfg)]),
    ),
)

visualizer = Visualizer(
    vis_backends=SwanlabVisBackend(init_kwargs={})
)  # init args can be found in https://docs.swanlab.cn/zh/guide_cloud/integration/integration-mmengine.html

runner = Runner(
    model=MMResNet50(),
    work_dir="./work_dir",
    train_dataloader=train_dataloader,
    optim_wrapper=dict(optimizer=dict(type=SGD, lr=0.001, momentum=0.9)),
    train_cfg=dict(by_epoch=True, max_epochs=5, val_interval=1),
    val_dataloader=val_dataloader,
    val_cfg=dict(),
    val_evaluator=dict(type=Accuracy),
    visualizer=visualizer,
)
runner.train()
