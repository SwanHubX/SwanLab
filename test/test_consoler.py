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

print("test myconsoler")
print("nihao swanlab")
