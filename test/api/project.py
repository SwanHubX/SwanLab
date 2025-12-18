"""
@author: Zhou Qiyang
@file: project.py
@time: 2025/12/17 10:45
@description: 用于测试api登录功能
"""
from typing import List
from unittest.mock import patch, MagicMock

import swanlab
from swanlab.api.model import Project

fake_projects: List[Project] = [
    Project({
        "cuid": "c1",
        "name": "proj-1",
        "path": "user/proj-1",
        "url": "https://example.com/@user/proj-1",
        "description": "desc-1",
        "visibility": "PUBLIC",
        "createdAt": "2025-01-01T00:00:00Z",
        "updatedAt": "2025-01-01T00:00:00Z",
        "projectLabels": [],
        "group": {"username": "user", "status": "ENABLED", "type": "TEAM"},
        "_count": {},
    })
]


def test_projects_uses_correct_params():
    with patch("swanlab.OpenApi") as MockOpenApi:
        mock_api = MagicMock()
        mock_api.projects.return_value = fake_projects
        MockOpenApi.return_value = mock_api

        # 这里导入或调用你的业务函数/脚本，它里面会执行 swanlab.OpenApi()
        api = swanlab.OpenApi()
        projects = api.projects(entity='bainiantest', sort=['create'], search="proj-1", detail=True)

        # 断言：OpenApi 被调用了一次
        MockOpenApi.assert_called_once()

        # 使用 mock 的 assert 方法判断 projects 是否为 fake_projects
        mock_api.projects.assert_called_once_with(entity='bainiantest', sort=['create'], search="proj-1", detail=True)
        assert projects == fake_projects
