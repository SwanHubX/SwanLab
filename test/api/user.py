"""
@author: Zhou QiYang
@file: user.py
@time: 2026/1/2 21:30
@description: OpenApi 用户相关测试代码
"""

import swanlab


def test_api_key():
    api = swanlab.Api()
    user = api.user

    print(user.__dict__)
    new_key = user.generate_api_key()
    print(new_key)
    print(user.api_keys)
    user.delete_api_key(api_key=new_key)
    print(user.api_keys)
