#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/8/22 14:52
@File: folder.py
@IDE: pycharm
@Description:
    文件夹上传模型
"""
from typing import List, Tuple
import click
from ..model import LaunchParser
from swanlab.error import ApiError
from swanlab.cli.utils import login_init_sid, UseTaskHttp, CosUploader, UploadBytesIO
import zipfile
from rich.progress import (
    BarColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
)
from swankit.log import FONT
from rich import print as rprint
from rich.filesize import decimal
from rich.text import Text
import time
import glob
import os
import io


class FolderParser(LaunchParser):
    """
    承担了上传文件夹的一系列任务
    """

    def __init__(self, config: dict, path: str):
        super().__init__(config, path)
        self.metadata = {"name": '', "desc": '', "combo": ''}
        self.spec = {"python": '', "entry": '', "volumes": [], "exclude": []}
        self.key = None
        """
        上传到cos的路径
        """
        self.api_key = None
        """
        用户的api_key
        """

    @classmethod
    def __type__(cls):
        return 'Folder'

    def __dict__(self) -> dict:
        data = {
            "src": self.key,
            "index": self.spec['entry'],
            "python": self.spec['python'],
            "conf": {"key": self.api_key},
            "name": self.metadata['name'],
        }
        if self.metadata.get("desc"):
            data["desc"] = self.metadata['desc']
        if self.metadata.get("combo"):
            data["combo"] = self.metadata['combo']
        if len(self.spec["volumes"]) > 0:
            data['datasets'] = [v['id'] for v in self.spec['volumes']]
        return data

    def parse_metadata(self, metadata: dict):
        name = self.should_be('metadata.name', metadata.get('name'), str)
        desc = self.should_be('metadata.desc', metadata.get('desc'), str, none=True)
        combo = self.should_be('metadata.combo', metadata.get('combo'), str, none=True)
        self.metadata['name'] = name
        self.metadata['desc'] = desc
        self.metadata['combo'] = combo

    def parse_spec(self, spec: dict):
        python = self.should_be('spec.python', spec.get('python'), str, none=True) or '3.10'
        python = self.should_in_values('spec.python', python, ['3.11', '3.10', '3.9', '3.8'])
        entry = self.should_be('spec.entry', spec.get('entry'), str, none=True) or 'train.py'
        self.should_file_exist('spec.entry', entry)
        volumes = self.should_be('spec.volumes', spec.get('volumes'), list, none=True) or []
        # NOTE 当前后端只支持一个volume
        import click

        if len(volumes) > 1:
            raise click.BadParameter('Only one volume is supported')

        for volume in volumes:
            self.should_equal_keys('volume', volume, ['name', 'id'])
            self.should_be('volume.name', volume.get('name'), str, none=True)
            self.should_be('volume.id', volume.get('id'), str)
        exclude = self.should_be('spec.exclude', spec.get('exclude'), list, none=True) or []
        [self.should_be('exclude', e, str) for e in exclude]

        self.spec['python'] = 'python' + python
        self.spec['entry'] = entry
        self.spec['volumes'] = volumes
        self.spec['exclude'] = exclude

    def walk(self, path: str = None) -> Tuple[List[str], List[str]]:
        """
        遍历path，生成文件列表，注意排除exclude中的文件
        此函数为递归调用函数
        返回所有命中的文件列表和排除的文件列表
        """
        path = path or self.dirpath
        all_files = glob.glob(os.path.join(path, '**'))
        exclude_files = []
        for g in self.spec['exclude']:
            efs = glob.glob(os.path.join(path, g))
            exclude_files.extend(efs)
        exclude_files = list(set(exclude_files))
        files = []
        for f in all_files:
            if os.path.isdir(f):
                if f in exclude_files:
                    continue
                fs, efs = self.walk(f)
                files.extend(fs)
                exclude_files.extend(efs)
            else:
                if f in exclude_files:
                    continue
                files.append(f)
        return files, exclude_files

    def zip(self, files: List[str]) -> io.BytesIO:
        """
        将walk得到的文件列表压缩到memory_file中
        """
        memory_file = io.BytesIO()
        progress = Progress(
            TextColumn("{task.description}", justify="left"),
            BarColumn(),
            "[progress.percentage]{task.percentage:>3.1f}%",
            "•",
            TimeRemainingColumn(),
        )
        z = zipfile.ZipFile(memory_file, "w", zipfile.ZIP_DEFLATED)
        with progress:
            for i in progress.track(range(len(files)), description=FONT.swanlab("Packing...  ")):
                arcname = os.path.relpath(files[i], start=self.dirpath)
                z.write(files[i], arcname)
        memory_file.seek(0)
        return memory_file

    def upload(self, memory_file: io.BytesIO):
        """
        上传压缩文件
        """
        val = memory_file.getvalue()
        client, sts = CosUploader.create()
        self.key = sts['prefix'] + "/tasks/" + f"{int(time.time() * 1000)}.zip"
        with UploadBytesIO(FONT.swanlab("Uploading..."), val) as buffer:
            client.upload_file_from_buffer(
                Bucket=sts['bucket'],
                Key=self.key,
                Body=buffer,
                MAXThread=5,
                MaxBufferSize=5,
                PartSize=1,
            )

    def run(self):
        # 剔除、压缩、上传、发布任务
        files, _ = self.walk()
        if len(files) == 0:
            raise click.BadParameter(self.dirpath + " is empty")
        login_info = login_init_sid()
        print(FONT.swanlab("Login successfully. Hi, " + FONT.bold(FONT.default(login_info.username))) + "!")
        self.api_key = login_info.api_key
        memory_file = self.zip(files)
        self.upload(memory_file)
        with UseTaskHttp() as http:
            try:
                http.post("/task", data=self.__dict__())
            except ApiError as e:
                if e.resp.status_code not in [404, 401]:
                    raise e
                elif e.resp.status_code == 404:
                    raise click.BadParameter("The dataset does not exist")
                else:
                    raise click.BadParameter("The combo does not exist")
        print(
            FONT.swanlab(
                f"Launch task successfully, use {FONT.bold(FONT.default('swanlab task list'))} to view the task"
            )
        )

    def dry_run(self):
        # 剔除、显示即将发布的任务的相关信息
        # 1. 任务名称
        # 2. 任务描述
        # 3. 任务套餐
        # 4. 任务python版本
        # 5. 任务入口文件路径
        # 6. 上传的任务文件夹路径
        # 7. 上传的任务文件夹中忽略的文件列表
        # 8. 绑定的数据卷信息
        _, exclude_files = self.walk()
        print(FONT.swanlab("This task will be launched:"))
        rprint("[bold]Name: [/bold]" + self.metadata['name'])
        rprint("[bold]Description: [/bold]" + self.metadata['desc'])
        rprint("[bold]Combo: [/bold]" + self.metadata['combo'] or 'Default')
        rprint("[bold]Python: [/bold]" + self.spec['python'])
        rprint("[bold]Entry: [/bold]" + self.spec['entry'])
        rprint("[bold]Folder: [/bold]" + self.dirpath)
        rprint("[bold]Excluded files: [/bold]")
        files = []
        for ef in exclude_files:
            if os.path.isdir(ef):
                suffix = '    📁 '
                index = 1  # 文件夹在前
                size = 0
                for path, dirs, ef_files in os.walk(ef):
                    for _ in ef_files:
                        fp = os.path.join(path, _)
                        size += os.path.getsize(fp)
            else:
                suffix = "    🐍 " if ef.endswith('.py') else "    📄 "
                index = 2
                size = os.path.getsize(ef)
            ef = os.path.relpath(ef, start=self.dirpath)
            files.append({'str': f'{ef} ({decimal(size)})', 'index': index, 'icon': Text(suffix)})
        files.sort(key=lambda x: x['str'])
        files.sort(key=lambda x: x['index'])
        for f in files:
            rprint(f['icon'] + f['str'])
        rprint("[bold]Volumes: [/bold]" + (str(self.spec['volumes']) if len(self.spec['volumes']) > 0 else 'None'))
