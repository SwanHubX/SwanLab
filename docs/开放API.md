# 开放API

> SwanLab SDK通过**开放API**提供对SwanLab云端服务的访问能力, 允许用户与SwanLab云端进行交互
> 
> 例如管理实验或项目, 获取个人资料等

## 用法

要使用开放API, 只需导入`swanlab.openapi`模块

```python
from swanlab import openapi as api

print(api.list_workspaces())
```

## 认证

开放API的认证与`swanlab.init()`的认证逻辑一致, 即在终端登录过, 或使用`swanlab.login()`登录
