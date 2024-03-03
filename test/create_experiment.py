import swanlab
import random

offset = random.random() / 5

run = swanlab.init(
    experiment_name="Example",
    description="这是一个机器学习模拟实验",
    config={
        "learning_rate": 0.01,
        "epochs": 20,
    },
)

# 模拟机器学习训练过程
for epoch in range(2, run.config.epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    swanlab.log({"loss": loss, "accuracy": acc})
