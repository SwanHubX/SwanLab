import swanlab
import time

swanlab.init()

for i in range(10):
    swanlab.log({"loss": i})
    time.sleep(1)

# the final chart
swanlab.log({"acc": 0.99})
