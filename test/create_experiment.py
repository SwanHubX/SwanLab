import swanlab
import time
import random
import numpy as np

epochs = 50
lr = 0.01
offset = random.random() / 5

swanlab.init(
    log_level="debug",
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
        arr1 = generate_random_nx3(1000)
        arr2 = generate_random_nx4(1000)
        arr3 = generate_random_nx6(1000)
        # swanlab.log(
        #     {
        #         "test/object3d1":
        #     },
        #     step=epoch,
        # )
        swanlab.log(
            {
                "test-object3d1": swanlab.Object3D(arr1, caption="3D点云数据"),
                "test-object3d2": swanlab.Object3D(arr2, caption="3D点云数据"),
                "test-object3d3": swanlab.Object3D(arr3, caption="3D点云数据"),
            },
            step=epoch,
        )
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    loss2 = 3**-epoch + random.random() / epoch + offset * 3
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    swanlab.log({"t/accuracy": acc, "loss": loss, "loss2": loss2})
    time.sleep(0.2)
