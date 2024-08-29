#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/8/22 14:56
@File: test_cli_launch.py
@IDE: pycharm
@Description:
    测试启动
    主要是配置文件的解析
"""
from swanlab.cli.commands import launcher as L
import tutils as T
import pytest
import click
import os


class TestParse:
    """
    测试解析配置文件，对应parse函数
    """

    def test_error_api_version(self):
        """
        测试错误的apiVersion
        """
        config = {
            'apiVersion': '12345',
        }
        with pytest.raises(click.UsageError) as e:
            L.parse(config, '')
        assert str(e.value) == 'Unknown api version: 12345'

    def test_error_kind(self):
        """
        测试错误的kind
        """
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'folder',  # 大小写敏感
        }
        with pytest.raises(click.UsageError) as e:
            L.parse(config, '')
        assert str(e.value) == 'Unknown kind: folder'

    def test_get_folder_v1(self):
        """
        测试成功
        """
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
        }
        parser = L.parse(config, '')
        assert isinstance(parser, L.parser.v1.FolderParser)


class TestFolderParserMetadata:
    """
    测试解析器解析metadata部分
    """

    def test_metadata_error_type(self):
        """
        测试错误的类型
        """
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            'metadata': 12345,
        }
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, '')
            parser.parse()
        assert str(e.value) == 'metadata should be <class \'dict\'>, not <class \'int\'>'

    def test_metadata_error_key(self):
        """
        测试错误的key
        """
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            'metadata': {
                'name': 'test',
                'error': 'test',
            },
        }
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, '')
            parser.parse()
        assert str(e.value) == 'Unknown key: metadata.error'

    def test_no_name(self):
        """
        测试没有name
        """
        config = {'apiVersion': 'swanlab/v1', 'kind': 'Folder', 'metadata': {}}
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, '')
            parser.parse()
        assert str(e.value) == 'metadata.name should not be None'

    def test_no_combo(self):
        """
        测试没有combo
        """
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            'metadata': {
                'name': 'test',
            },
        }
        parser = L.parse(config, '')
        parser.parse_metadata(config['metadata'])
        assert parser.metadata['combo'] == None  # noqa

    def test_no_desc(self):
        """
        测试没有desc
        """
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            'metadata': {
                'name': 'test',
                'combo': 'test',
            },
        }
        parser = L.parse(config, '')
        parser.parse_metadata(config['metadata'])
        assert parser.metadata['desc'] is None  # noqa
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            'metadata': {
                'name': 'test',
                'combo': 'test',
                'desc': 'test',
            },
        }
        parser = L.parse(config, '')
        parser.parse_metadata(config['metadata'])
        assert parser.metadata['desc'] == 'test'  # noqa


class TestFolderParserSpec:
    """
    测试解析器解析spec部分
    """

    def test_spec_error_type(self):
        """
        测试错误的类型
        """
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            "metadata": {"name": "test", "desc": "test", "combo": "test"},
            'spec': 12345,
        }
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, '')
            parser.parse()
        assert str(e.value) == 'spec should be <class \'dict\'>, not <class \'int\'>'

    def test_spec_error_key(self):
        """
        测试错误的key
        """
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            "metadata": {"name": "test", "desc": "test", "combo": "test"},
            'spec': {
                'name': 'test',
            },
        }
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, '')
            parser.parse()
        assert str(e.value) == 'Unknown key: spec.name'

    @staticmethod
    def mock_entry(name='train.py') -> str:
        """
        模拟一个文件
        """
        with open(os.path.join(T.TEMP_PATH, name), 'w') as f:
            f.write('print("hello")')
        return os.path.join(T.TEMP_PATH, 'swanlab.yaml')

    def test_python(self):
        f = self.mock_entry()
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            "metadata": {"name": "test", "desc": "test", "combo": "test"},
            'spec': {
                'entry': 'train.py',
                'python': '3.8',
            },
        }
        parser = L.parse(config, f)
        parser.parse()
        assert parser.spec['python'] == 'python3.8'  # noqa

    def test_error_python(self):
        f = self.mock_entry()
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            "metadata": {"name": "test", "desc": "test", "combo": "test"},
            'spec': {
                'entry': 'train.py',
                'python': '3.7',
            },
        }
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, f)
            parser.parse()
        assert str(e.value) == 'spec.python should be in [\'3.11\', \'3.10\', \'3.9\', \'3.8\'], not 3.7'

    def test_no_python(self):
        """
        测试没有python
        """
        f = self.mock_entry()
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            "metadata": {"name": "test", "desc": "test", "combo": "test"},
            'spec': {
                'entry': 'train.py',
            },
        }
        parser = L.parse(config, f)
        parser.parse()
        assert parser.spec['python'] == 'python3.10'  # noqa

    def test_no_entry(self):
        """
        测试没有entry
        """
        f = self.mock_entry()
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            "metadata": {"name": "test", "desc": "test", "combo": "test"},
            'spec': {
                'python': '3.8',
            },
        }
        parser = L.parse(config, f)
        parser.parse()
        parser.spec['entry'] == 'train.py'  # noqa

    def test_no_entry_file(self):
        """
        测试entry文件不存在
        """
        f = self.mock_entry('test.py')
        err_entry = '1234'
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            "metadata": {"name": "test", "desc": "test", "combo": "test"},
            'spec': {
                'entry': err_entry,
            },
        }
        with pytest.raises(click.FileError) as e:
            parser = L.parse(config, f)
            parser.parse()
        assert str(e.value) == f'spec.entry not found: {os.path.join(T.TEMP_PATH, err_entry)}'

    def test_entry_not_file(self):
        """
        测试entry不是文件
        """
        f = self.mock_entry()
        err_entry = 'test'
        os.mkdir(os.path.join(T.TEMP_PATH, err_entry))
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            "metadata": {"name": "test", "desc": "test", "combo": "test"},
            'spec': {
                'entry': err_entry,
            },
        }
        with pytest.raises(click.FileError) as e:
            parser = L.parse(config, f)
            parser.parse()
        assert str(e.value) == f'spec.entry should be a file: {os.path.join(T.TEMP_PATH, err_entry)}'

    def test_volumes(self):
        f = self.mock_entry()
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            "metadata": {"name": "test", "desc": "test", "combo": "test"},
            'spec': {
                'entry': 'train.py',
                'volumes': [
                    {'name': 'test', 'id': 'test'},
                ],
            },
        }
        parser = L.parse(config, f)
        parser.parse()
        assert parser.spec['volumes'] == [{'name': 'test', 'id': 'test'}]  # noqa
        config['spec']['volumes'].append({'name': 'test', 'id': '2'})
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, f)
            parser.parse()
        assert str(e.value) == 'Only one volume is supported'
        config['spec']['volumes'] = [{'id': 'test'}]  # noqa
        parser = L.parse(config, f)
        parser.parse()
        assert parser.spec['volumes'] == [{'id': 'test'}]  # noqa
        config['spec']['volumes'] = [{'name': 'test'}]  # noqa
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, f)
            parser.parse()
        assert str(e.value) == 'volume.id should not be None'

    def test_error_volumes_type(self):
        f = self.mock_entry()
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            "metadata": {"name": "test", "desc": "test", "combo": "test"},
            'spec': {
                'entry': 'train.py',
                'volumes': 'test',
            },
        }
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, f)
            parser.parse()
        assert str(e.value) == 'spec.volumes should be <class \'list\'>, not <class \'str\'>'
        config['spec']['volumes'] = [{'name': 'test', 'id': 123}]  # noqa
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, f)
            parser.parse()
        assert str(e.value) == 'volume.id should be <class \'str\'>, not <class \'int\'>'
        config['spec']['volumes'] = [{'name': 'test', 'id': 'test', 'error': 1}]  # noqa
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, f)
            parser.parse()
        assert str(e.value) == 'Unknown key: volume.error'

    def test_exclude(self):
        f = self.mock_entry()
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            "metadata": {"name": "test", "desc": "test", "combo": "test"},
            'spec': {
                'entry': 'train.py',
                'exclude': ['test'],
            },
        }
        parser = L.parse(config, f)
        parser.parse()
        assert parser.spec['exclude'] == ['test']  # noqa

    def test_exclude_error_type(self):
        f = self.mock_entry()
        config = {
            'apiVersion': 'swanlab/v1',
            'kind': 'Folder',
            "metadata": {"name": "test", "desc": "test", "combo": "test"},
            'spec': {
                'entry': 'train.py',
                'exclude': 'test',
            },
        }
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, f)
            parser.parse()
        assert str(e.value) == 'spec.exclude should be <class \'list\'>, not <class \'str\'>'
        config['spec']['exclude'] = ['test', 123]  # noqa
        with pytest.raises(click.BadParameter) as e:
            parser = L.parse(config, f)
            parser.parse()
        assert str(e.value) == 'exclude should be <class \'str\'>, not <class \'int\'>'
