import xgboost as xgb
from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report
import swanlab
from swanlab.integration.xgboost import SwanLabCallback

# 初始化swanlab
swanlab.init(project="xgboost-breast-cancer", config={
    "learning_rate": 0.1,
    "max_depth": 3,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "num_round": 100
})

# 加载数据集
data = load_breast_cancer()
X = data.data
y = data.target

# 将数据集分为训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 转换为DMatrix格式，这是XGBoost的内部数据格式
dtrain = xgb.DMatrix(X_train, label=y_train)
dtest = xgb.DMatrix(X_test, label=y_test)

# 设置参数
params = {
    'objective': 'binary:logistic',  # 二分类任务
    'max_depth': 3,                  # 树的最大深度
    'eta': 0.1,                      # 学习率
    'subsample': 0.8,                # 样本采样比例
    'colsample_bytree': 0.8,         # 特征采样比例
    'eval_metric': 'logloss'         # 评估指标
}

# 训练模型
num_round = 100  # 迭代次数
bst = xgb.train(params, dtrain, num_round, evals=[(dtrain, 'train'), (dtest, 'test')], callbacks=[SwanLabCallback()])

# 进行预测
y_pred = bst.predict(dtest)
y_pred_binary = [round(value) for value in y_pred]  # 将概率转换为二分类结果

# 评估模型
accuracy = accuracy_score(y_test, y_pred_binary)
print(f"Accuracy: {accuracy:.4f}")

# 打印分类报告
print("Classification Report:")
print(classification_report(y_test, y_pred_binary, target_names=data.target_names))

# 保存模型
bst.save_model('xgboost_model.model')

# 结束swanlab会话
swanlab.finish()