#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/8/22 14:52
@File: folder.py
@IDE: pycharm
@Description:
    文件夹上传模型
"""
from ..model import LaunchParser


class FolderParser(LaunchParser):
    """
    承担了上传文件夹的一系列任务
    """

    def __init__(self, config: dict, path: str):
        super().__init__(config, path)
        self.metadata = {
            "name": None,
            "desc": None,
            "combo": None
        }
        self.spec = {
            "python": None,
            "entry": None,
            "volumes": [],
            "exclude": []
        }

    @classmethod
    def __type__(cls):
        return 'Folder'

    def __dict__(self) -> dict:
        pass

    def parse_metadata(self, metadata: dict):
        name = self.should_be('metadata.name', metadata.get('name'), str)
        desc = self.should_be('metadata.desc', metadata.get('desc'), str, none=True)
        combo = self.should_be('metadata.combo', metadata.get('combo'), str)
        self.metadata['name'] = name
        self.metadata['desc'] = desc
        self.metadata['combo'] = combo

    def parse_spec(self, spec: dict):
        python = self.should_be('spec.python', spec.get('python'), str, none=True) or '3.10'
        python = self.should_in_values('spec.python', python, ['3.11', '3.10', '3.9', '3.8'])
        entry = self.should_be('spec.entry', spec.get('entry'), str) or 'train.py'
        self.should_file_exist('spec.entry', entry)
        volumes = self.should_be('spec.volumes', spec.get('volumes'), list, none=True) or []
        # NOTE 当前后端只支持一个volume
        import click
        if len(volumes) > 1:
            raise click.BadParameter('Only one volume is supported')

        for volume in volumes:
            self.should_equal_keys('volume', volume, ['name', 'id'])
            self.should_be('volume.name', volume.get('name'), str)
            self.should_be('volume.id', volume.get('id'), str)
        exclude = self.should_be('spec.exclude', spec.get('exclude'), list, none=True) or []
        [self.should_be('exclude', e, str) for e in exclude]

        self.spec['python'] = python
        self.spec['entry'] = entry
        self.spec['volumes'] = volumes
        self.spec['exclude'] = exclude

    def run(self):
        pass

    def dry_run(self):
        pass
