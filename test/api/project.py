"""
@author: Zhou Qiyang
@file: project.py
@time: 2025/12/17 10:45
@description: 用于测试api登录功能
"""

from unittest.mock import patch

import swanlab

# 测试数据
fake_projects_raw = [
    {
        "total": 10,
        "pages": 1,
        "size": 10,
        "list": [
            {
                "cuid": "c1",
                "name": "proj-1",
                "path": "user/proj-1",
                "url": "https://example.com/@user/proj-1",
                "description": "desc-1",
                "visibility": "PUBLIC",
                "createdAt": "2025-01-01T00:00:00Z",
                "updatedAt": "2025-01-01T00:00:00Z",
                "projectLabels": [{"name": "Nvidia"}],
                "group": {"username": "user", "status": "ENABLED", "type": "TEAM"},
                "_count": {},
            }
        ],
    }
]


def test_projects_uses_correct_params():
    # patch: client.get 返回 fake_projects_raw
    with patch("swanlab.core_python.client.Client.get", return_value=fake_projects_raw) as mock_get:
        api = swanlab.OpenApi()

        result = api.projects(entity="bainiantest", detail=True)

        # 断言返回值
        # 1. 字符串类型的字段
        raw_list = fake_projects_raw[0]["list"]
        fields = [
            "cuid",
            "name",
            "path",
            "url",
            "description",
            "visibility",
            "createdAt",
            "updatedAt",
        ]
        for field in fields:
            assert [getattr(p, field) for p in result] == [r[field] for r in raw_list]

        # 2. labels
        assert [[l.name for l in p.projectLabels] for p in result] == [
            [l["name"] for l in r["projectLabels"]] for r in raw_list
        ]

        # 3. group
        for field in ["username", "type", "status"]:
            assert [getattr(p.group, field) for p in result] == [r["group"][field] for r in raw_list]

        # 4. count
        assert [p.count for p in result] == [r["_count"] for r in raw_list]
