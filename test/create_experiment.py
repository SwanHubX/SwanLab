import swanlab
import time

swanlab.init(log_level="debug")

for i in range(100):
    swanlab.log({"lllloss": i})
    # time.sleep(0.1)

for i in range(300):
    swanlab.log({"accccc": "xfds"})
    # time.sleep(0.1)
