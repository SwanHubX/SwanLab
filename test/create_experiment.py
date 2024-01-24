import swanlab
import time
import random

swanlab.init(log_level="debug")

for i in range(100):
    swanlab.log({"loss": random.choice([0.1, 0.2, 0.3, 0.4, "nan"]), "acc": i})
    time.sleep(0.1)

for i in range(300):
    swanlab.log({"acc": 1})
    time.sleep(0.1)
