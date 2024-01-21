import swanlab
import time

swanlab.init()

for i in range(10):
    swanlab.log({"loss": i})
    time.sleep(0.5)

for i in range(30):
    swanlab.log({"acc": 1})
    time.sleep(0.5)
