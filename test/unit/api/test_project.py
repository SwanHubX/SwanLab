from unittest.mock import patch, MagicMock

from swanlab.api.projects import Projects
from swanlab.core_python import Client
from swanlab.package import get_host_web
from utils import create_project_data


def test_projects():
    """测试能否分页获取所有项目"""
    with patch('swanlab.api.projects.get_workspace_projects') as mock_get_projects:
        total_pages = 4
        page_size = 20

        def side_effect(*args, **kwargs):
            return create_project_data(size=page_size, pages=kwargs.get("page", 1), total=page_size * total_pages)

        mock_get_projects.side_effect = side_effect

        mock_projects = Projects(
            MagicMock(spec=Client),
            web_host=get_host_web(),
            workspace='test_user',
        )
        projects = list(mock_projects)
        assert len(projects) == page_size * total_pages
        assert mock_get_projects.call_count == total_pages
