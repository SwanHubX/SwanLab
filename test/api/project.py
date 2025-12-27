"""
@author: Zhou Qiyang
@file: project.py
@time: 2025/12/17 10:45
@description: 用于测试api登录功能
"""

from unittest.mock import patch

import swanlab


def example_code():
    api = swanlab.OpenApi()
    projects = api.projects(workspace='workspace')
    for project in projects:
        print(project.__dict__)


# 测试数据
def make_fake_projects(start, count):
    return [
        {
            "cuid": f"c{n}",
            "name": f"proj-{n}",
            "path": f"user/proj-{n}",
            "url": f"https://dev001.swanlab.cn/@user/proj-{n}",
            "description": f"desc-{n}",
            "visibility": "PUBLIC",
            "createdAt": "2025-01-01T00:00:00Z",
            "updatedAt": "2025-01-01T00:00:00Z",
            "projectLabels": [{"name": "Nvidia"}],
            "group": {"username": "user", "status": "ENABLED", "type": "TEAM"},
            "_count": {},
        }
        for n in range(start, start + count)
    ]


performance_test_projects = [
    [
        {
            "total": 40,
            "pages": 1,
            "size": 20,
            "list": make_fake_projects(0, 20),
        },
    ],
    [
        {
            "total": 40,
            "pages": 2,
            "size": 20,
            "list": make_fake_projects(20, 20),
        },
    ],
]

params_test_projects = [
    {
        "total": 20,
        "pages": 1,
        "size": 20,
        "list": make_fake_projects(0, 20),
    },
]


# 性能测试：是否按照当前遍历的项目动态获取
def test_api_projects_performance():
    # patch: client.get 返回 fake_projects_raw
    with patch("swanlab.core_python.client.Client.get", side_effect=performance_test_projects) as mock_get:
        api = swanlab.OpenApi()

        result = api.projects(workspace="user", detail=True)

        # 断言请求调用次数
        assert mock_get.call_count == 0
        for project in result:
            if project.name == "proj-19":
                assert mock_get.call_count == 1
            if project.name == "proj-20":
                assert mock_get.call_count == 2


# 功能测试：获取到的项目的属性是否齐全且正确
def test_api_projects_params():
    with patch("swanlab.core_python.client.Client.get", return_value=params_test_projects):
        api = swanlab.OpenApi()

        result = api.projects(workspace="user", detail=True)

        # 1. 字符串类型的字段
        raw_list = params_test_projects[0]["list"]
        fields = {
            "name": "name",
            "path": "path",
            "url": "url",
            "description": "description",
            "visibility": "visibility",
            "created_at": "createdAt",
            "updated_at": "updatedAt",
        }
        for field in fields:
            assert [getattr(p, field) for p in result] == [r[fields[field]] for r in raw_list]

        # 1.1 workspace
        assert [p.workspace for p in result] == [r["group"]["username"] for r in raw_list]

        # 2. labels
        assert [[l.__str__() for l in p.labels] for p in result] == [
            [l["name"] for l in r["projectLabels"]] for r in raw_list
        ]
        assert [[l.name for l in p.labels] for p in result] == [
            [l["name"] for l in r["projectLabels"]] for r in raw_list
        ]

        # 3. count
        assert [p.count for p in result] == [r["_count"] for r in raw_list]
