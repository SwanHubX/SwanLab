from swanlab.log import init_consoler
from swanlab.env import swc
import swanlab as sw

# 迭代次数
epochs = 200
# 学习率
lr = 0.01
# 随机偏移量
offset = 0.5

# 创建一个实验
sw.init(
    description="this is a test experiment",
    config={
        "learning_rate": lr,
        "epochs": epochs,
    },
)

swc.init(swc.getcwd(), "train")

init_consoler(swc.console_folder)

print("test myconsoler")
print("nihao swanlab")
