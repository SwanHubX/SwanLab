import swanlab
import time

swanlab.init(
    log_level="debug",
    config={
        "test": 1,
        "debug": 1,
        "verbose": 1,
    },
)

for i in range(100):
    swanlab.log({"loss": i})
    time.sleep(0.1)

for i in range(300):
    swanlab.log({"acc": 1})
    time.sleep(0.1)
