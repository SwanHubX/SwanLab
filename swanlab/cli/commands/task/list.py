#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/19 14:09
@File: status.py
@IDE: pycharm
@Description:
    列出任务状态
"""
import time

import click
from typing import List
from .utils import login_init_sid
from swankit.log import FONT
from rich.layout import Layout
from datetime import datetime
from rich.panel import Panel
from rich.table import Table
from rich.live import Live
from swanlab.package import get_experiment_url


@click.command()
@click.option(
    "--max-num",
    "-m",
    default=10,
    nargs=1,
    type=click.IntRange(1, 100),
    help="The maximum number of tasks to display, default by 10, maximum by 100",
)
def list(max_num: int):  # noqa
    # 获取访问凭证，生成http会话对象
    login_info = login_init_sid()
    print(FONT.swanlab("You are logging as " + FONT.bold(FONT.default(login_info.username))) + ".")
    # 获取任务列表
    ltm = ListTasksModel(num=max_num, username=login_info.username)
    layout = ListTaskLayout(ltm)
    layout.show()


class ListTasksModel:
    class TaskListModel:
        """
        获取到的任务列表模型
        """

        def __init__(self, username: str, task: dict, ):
            self.username = username
            self.name = task["name"]
            """
            任务名称
            """
            self.python = task["python"]
            """
            任务的python版本
            """
            self.project_name = task.get("pName", None)
            """
            项目名称
            """
            self.experiment_id = task.get("eId", None)
            """
            实验ID
            """
            self.created_at = task["createdAt"]
            self.started_at = task.get("startedAt", None)
            self.finished_at = task.get("finishedAt", None)
            self.status = task["status"]
            self.msg = task.get("msg", None)

        @property
        def url(self):
            if self.project_name is None or self.experiment_id is None:
                return None
            return get_experiment_url(self.username, self.project_name, self.experiment_id)

    def __init__(self, num: int, username: str):
        """
        :param num: 最大显示的任务数
        """
        self.num = num
        self.username = username

    def __dict__(self):
        return {"num": self.num}

    def list(self) -> List[TaskListModel]:
        # TODO 部署完毕接入http
        import requests
        resp = requests.get(
            "http://172.16.42.24:1323/api/task",
            json=self.__dict__(),
            headers={"payload": '{"uid": 1, "username": "' + str(self.username) + '"}'}
        )
        if resp.status_code != 200:
            raise ValueError(f"Error: {resp.json()}")
        return [self.TaskListModel(self.username, task) for task in resp.json()]

    def table(self):
        st = Table(
            expand=True,
            show_header=True,
            header_style="bold",
            title="[magenta][b]Now Tasks![/b]",
            highlight=True,
        )
        st.add_column("Task Name", justify="right")
        st.add_column("Status", justify="center")
        st.add_column("URL", justify="center")
        st.add_column("Python Version", justify="center"),
        st.add_column("Created Time", justify="center")
        st.add_column("Started Time", justify="center")
        st.add_column("Finished Time", justify="center")
        for tlm in self.list():
            st.add_row(
                tlm.name,
                tlm.status,
                tlm.url,
                tlm.python,
                tlm.created_at,
                tlm.started_at,
                tlm.finished_at,
            )
        time.sleep(10)
        return st


class ListTaskHeader:
    """
    Display header with clock.
    """

    @staticmethod
    def __rich__() -> Panel:
        grid = Table.grid(expand=True)
        grid.add_column(justify="center", ratio=1)
        grid.add_column(justify="right")
        grid.add_row(
            "[b]SwanLab[/b] task dashboard",
            datetime.now().ctime().replace(":", "[blink]:[/]"),
        )
        return Panel(grid, style="red on black")


class ListTaskLayout:
    """
    任务列表展示
    """

    def __init__(self, ltm: ListTasksModel):
        self.layout = Layout()
        self.layout.split(
            Layout(name="header", size=3),
            Layout(name="main")
        )
        self.layout["header"].update(ListTaskHeader())
        self.layout["main"].update(Panel(ltm.table(), border_style="magenta"))
        self.ltm = ltm

    def show(self):
        with Live(self.layout, refresh_per_second=10, screen=True) as live:
            now = time.time()
            while True:
                time.sleep(1)
                self.layout["header"].update(ListTaskHeader())
                if time.time() - now > 5:
                    now = time.time()
                    self.layout["main"].update(Panel(self.ltm.table(), border_style="magenta"))
                live.refresh()
