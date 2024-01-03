from test1 import test1
from swanlab.log import swanlog as sl
import swanlab as sw

# 迭代次数
epochs = 200
# 学习率
lr = 0.01

# 创建一个实验
sw.init(
    description="this is a test experiment",
    config={
        "learning_rate": lr,
        "epochs": epochs,
    },
)
sl.init("output.log", "debug")


print(sl.isRunning)
sl.setCollectionLevel("info")
print(sl.getCollectionLevel())
sl.debug("Watch out!")
sl.info("I told you so")
sl.warning("I told you so")
sl.error("I told you so")
sl.critical("I told you so")
# test1()

# sl.setSuccess()
