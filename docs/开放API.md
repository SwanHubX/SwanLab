# 开放API

> SwanLab SDK通过**开放API**提供对SwanLab云端服务的访问能力, 允许用户与SwanLab云端进行交互
> 
> 例如管理实验或项目, 获取个人资料等

## 用法

要使用开放API, 只需导入`swanlab.OpenApi`模块

```python
from swanlab import OpenApi

my_api = OpenApi() # 使用之前的登录信息
print(my_api.list_workspaces())

other_api = OpenApi(api_key='other_api_key') # 使用另一个账户的key
print(other_api.list_workspaces())
```

## 认证

- 如果显式提供`api_key`, 则使用该`api_key`进行认证
- 否则开放API的认证与`swanlab.init()`的认证逻辑一致, 即
  - 曾在终端登录过
  - 显式使用`swanlab.login()`登录
