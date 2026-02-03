import csv
from io import StringIO
from typing import Dict, List

import nanoid

from swanlab.core_python.api.type import RunType, ProjResponseType, ProjectType
from swanlab.package import get_host_web


def create_run_type_data(cuid=None) -> RunType:
    """
    创建 RunType 类型的测试数据
    """
    return {
        'cuid': cuid if cuid is not None else nanoid.generate('0123456789', 10),
        'name': '',
        'createdAt': '',
        'description': '',
        'labels': [],
        'profile': {'config': {}, 'scalar': {'loss': [], 'accuracy': []}},
        'state': 'FINISHED',
        'cluster': '',
        'job': '',
        'runtime': '',
        'user': {'username': '', 'name': ''},
        'rootExpId': None,
        'rootProId': None,
    }


def create_nested_exps(groups: int = 2, num_per_group: int = 2) -> Dict:
    """
    创建嵌套的实验数据（用于模拟 get_project_experiments 的返回值）

    :param groups: 分组数量
    :param num_per_group: 每组（页）中的实验数量
    :return: 分组情况下 Dict 格式的数据
    """
    result = {}
    for i in range(groups):
        group_key = f'group_{i}'
        runs = []
        for j in range(num_per_group):
            run = create_run_type_data(cuid=f'exp_{i}_{j}')
            runs.append(run)
        result[group_key] = runs
    return result


def create_project_data(page: int = 1, total: int = 20) -> ProjResponseType:
    """
    创建分页项目数据（用于模拟 get_workspace_projects 的返回值）

    :param page: 当前页数
    :param total: 项目总数
    :return: ProjResponseType 格式的数据
    """
    page_size = 20
    pages = (total + page_size - 1) // page_size
    project_list: List[ProjectType] = []

    for j in range(page_size):
        project: ProjectType = {
            'cuid': f'proj_{page}_{j}',
            'name': f'project_{page}_{j}',
            'path': f'test_user/project_{page}_{j}',
            'url': f'{get_host_web()}/test_user/project_{page}_{j}',
            'description': '',
            'visibility': 'PRIVATE',
            'createdAt': '',
            'updatedAt': '',
            'group': {'workspace': 'test_user'},
            'projectLabels': [],
            '_count': {},
        }
        project_list.append(project)

    return {
        'list': project_list,
        'size': page_size,
        'pages': pages,
        'total': total,
    }


def create_user_data(page: int = 1, total: int = 20) -> Dict:
    """
    创建分页用户数据（用于模拟 get_users 的返回值）

    :param page: 当前页数
    :param total: 用户总数
    :return: 包含 list, pages, total 等字段的字典
    """
    page_size = 20
    pages = (total + page_size - 1) // page_size
    user_list = []

    for j in range(page_size):
        user_list.append({
            'username': f'user_{ (page - 1) * page_size + j }'
        })

    return {
        'list': user_list,
        'pages': pages,
        'total': total,
    }


def create_csv_data(step_values, metric_name, metric_values):
    """
    创建 CSV 格式的数据
    :param step_values: step 列的值列表
    :param metric_name: 指标名称
    :param metric_values: 指标值列表
    :return: CSV 格式的字节数据
    """
    output = StringIO()
    writer = csv.writer(output)
    writer.writerow(['step', metric_name])
    for step, value in zip(step_values, metric_values):
        writer.writerow([step, value])
    csv_text = output.getvalue()
    return csv_text.encode('utf-8')
