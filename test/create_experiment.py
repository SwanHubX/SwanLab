import swanlab
import time

swanlab.init(
    log_level="info",
    config={
        "test": 1,
        "debug": "11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111",
        "verbose": 1,
    },
    logggings=True,
)

for i in range(30):
    swanlab.log({"acc": i, "hhh/loss": "123"})
    time.sleep(0.5)
