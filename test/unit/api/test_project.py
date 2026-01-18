from unittest.mock import patch, MagicMock
import responses
from responses import registries

from swanlab.api.projects import Projects
from swanlab.core_python import Client
from swanlab.package import get_host_web, get_host_api
from utils import create_project_data


def test_projects():
    """测试能否分页获取所有项目"""
    with patch('swanlab.api.projects.get_workspace_projects') as mock_get_projects:
        total_pages = 4
        page_size = 20
        mock_get_projects.return_value = create_project_data(
            pages=mock_get_projects.call_count + 1, size=page_size, total=page_size * total_pages
        )

        mock_projects = Projects(
            MagicMock(spec=Client),
            web_host=get_host_web(),
            workspace='test_user',
        )
        projects = list(mock_projects)
        assert len(projects) == page_size * total_pages
        assert mock_get_projects.call_count == total_pages
