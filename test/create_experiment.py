import swanlab
import random
import numpy as np

epochs = 50
lr = 0.01
offset = random.random() / 5

run = swanlab.init(
    experiment_name="Example",
    description="这是一个机器学习模拟实验",
    config={
        "epochs": epochs,
        "learning_rate": lr,
        "test": 1,
        "debug": "这是一串" + "很长" * 100 + "的字符串",
        "verbose": 1,
    },
    logggings=True,
)


def generate_random_nx3(n):
    """生成形状为nx3的随机数组"""
    return np.random.rand(n, 3)


def generate_random_nx4(n):
    """生成形状为nx4的随机数组，最后一列是[1,14]范围内的整数分类"""
    xyz = np.random.rand(n, 3)
    c = np.random.randint(1, 15, size=(n, 1))
    return np.hstack((xyz, c))


def generate_random_nx6(n):
    """生成形状为nx6的随机数组，包含RGB颜色"""
    xyz = np.random.rand(n, 3)
    rgb = np.random.rand(n, 3)  # RGB颜色值也可以是[0,1]之间的随机数
    rgb = (rgb * 255).astype(np.uint8)  # 转换为[0,255]之间的整数
    return np.hstack((xyz, rgb))


for epoch in range(2, epochs):
    if epoch % 10 == 0:

        # swanlab.log(
        #     {
        #         "test/object3d1":
        #     },
        #     step=epoch,
        # )
        swanlab.log(
            {
                "test-object3d1": swanlab.Object3D("./assets/bunny.obj", caption="bunny-obj"),
                "test-object3d2": swanlab.Object3D("./assets/test1.pts.json", caption="test1-pts"),
            },
            step=epoch,
        )
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    swanlab.log({"loss": loss, "accuracy": acc})
