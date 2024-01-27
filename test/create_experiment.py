import swanlab
import time
import random

epochs = 100
lr = 0.01
offset = random.random() / 5

swanlab.init(
    log_level="debug",
    config={
        "epochs": epochs,
        "learning_rate": lr,
        "test": 1,
        "debug": "11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111",
        "verbose": 1,
    },
    logggings=True,
)

for epoch in range(2, epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    swanlab.log({"accuracy": acc, "loss": loss})
    time.sleep(1)
    time.sleep(0.5)
