from unittest.mock import patch, MagicMock

from swanlab.api.experiments import Experiments
from swanlab.core_python import Client
from tutils.setup import mock_login_info
from utils import create_nested_exps


def test_folded_exps():
    """测试嵌套的实验数据能够正确展平"""
    mock_exps = Experiments(
        MagicMock(spec=Client),
        path='test_user/test-project',
        login_info=mock_login_info(),
    )

    nested_data = create_nested_exps(groups=2, num_per_group=2)
    with patch('swanlab.api.experiments.get_project_experiments') as mock_get_exps:
        mock_get_exps.return_value = nested_data
        experiments = list(mock_exps)
        assert len(experiments) == 4
