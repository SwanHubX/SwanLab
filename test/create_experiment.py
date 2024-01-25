import swanlab
import time

swanlab.init(
    log_level="info",
    config={
        "test": 1,
        "debug": 1,
        "verbose": 1,
    },
    logggings=True,
)

for i in range(10):
    swanlab.log({"acc": i, "loss": "123"})
