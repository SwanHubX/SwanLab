import lightgbm as lgb
from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import swanlab
from swanlab.integration.lightgbm import SwanLabCallback

# Step 1: Initialize swanlab
swanlab.init(project="lightgbm-example", name="breast-cancer-classification")

# Step 2: Load the dataset
data = load_breast_cancer()
X = data.data
y = data.target

# Step 3: Split the data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Step 4: Create LightGBM datasets
train_data = lgb.Dataset(X_train, label=y_train)
test_data = lgb.Dataset(X_test, label=y_test, reference=train_data)

# Step 5: Set parameters
params = {
    'objective': 'binary',
    'metric': 'binary_logloss',
    'boosting_type': 'gbdt',
    'num_leaves': 31,
    'learning_rate': 0.05,
    'feature_fraction': 0.9
}

# Step 6: Train the model with swanlab callback
num_round = 100
gbm = lgb.train(
    params,
    train_data,
    num_round,
    valid_sets=[test_data],
    callbacks=[SwanLabCallback()]
)

# Step 8: Make predictions
y_pred = gbm.predict(X_test, num_iteration=gbm.best_iteration)
y_pred_binary = [1 if p >= 0.5 else 0 for p in y_pred]

# Step 9: Evaluate the model
accuracy = accuracy_score(y_test, y_pred_binary)
print(f"模型准确率: {accuracy:.4f}")
swanlab.log({"accuracy": accuracy})

# Step 10: Save the model locally
gbm.save_model('lightgbm_model.txt')

# Step 11: Load the model and predict again
bst_loaded = lgb.Booster(model_file='lightgbm_model.txt')
y_pred_loaded = bst_loaded.predict(X_test)
y_pred_binary_loaded = [1 if p >= 0.5 else 0 for p in y_pred_loaded]

# Step 12: Evaluate the loaded model
accuracy_loaded = accuracy_score(y_test, y_pred_binary_loaded)
print(f"加载模型后的准确率: {accuracy_loaded:.4f}")
swanlab.log({"accuracy_loaded": accuracy_loaded})

# Step 13: Finish the swanlab run
swanlab.finish()