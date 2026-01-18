import csv
from io import StringIO
from typing import Dict, List
import nanoid
from swanlab.core_python.api.type import RunType, ProjResponseType, ProjectType
from swanlab.package import get_host_web


def create_run_type_data(**kwargs) -> RunType:
    """
    创建 RunType 类型的测试数据
    """
    return {
        'cuid': kwargs.get('cuid', nanoid.generate('0123456789', 10)),
        'name': kwargs.get('name', ''),
        'createdAt': kwargs.get('createdAt', ''),
        'description': kwargs.get('description', ''),
        'labels': kwargs.get('labels', []),
        'profile': kwargs.get('profile', {'config': {}, 'scalar': {'loss': [], 'accuracy': []}}),
        'state': kwargs.get('state', 'FINISHED'),
        'cluster': kwargs.get('cluster', ''),
        'job': kwargs.get('job', ''),
        'runtime': kwargs.get('runtime', ''),
        'user': kwargs.get('user', {'username': '', 'name': ''}),
        'rootExpId': kwargs.get('rootExpId', None),
        'rootProId': kwargs.get('rootProId', None),
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


def create_project_data(size: int = 20, pages: int = 1, total: int = 20) -> ProjResponseType:
    """
    创建单页项目数据（用于模拟 get_workspace_projects 的返回值，仅生成一页）

    :param size: 项目数
    :param pages: 当前页数
    :param total: 项目总数
    :return: ProjResponseType 格式的数据
    """
    project_list: List[ProjectType] = []

    for j in range(size):
        project: ProjectType = {
            'cuid': f'proj_0_{j}',
            'name': f'project_0_{j}',
            'path': f'test_user/project_0_{j}',
            'url': f'{get_host_web()}/test_user/project_0_{j}',
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
        'size': size,
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
    # 写入表头：step, metric_name
    writer.writerow(['step', metric_name])
    # 写入数据
    for step, value in zip(step_values, metric_values):
        writer.writerow([step, value])
    # 获取文本内容并编码为字节
    csv_text = output.getvalue()
    return csv_text.encode('utf-8')
