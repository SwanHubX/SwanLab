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
from rich.layout import Layout
from datetime import datetime
from rich.panel import Panel
from rich.table import Table
from rich.live import Live
from .utils import TaskModel
from swanlab.cli.utils import login_init_sid, UseTaskHttp


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
    """
    List tasks
    """
    # 获取访问凭证，生成http会话对象
    login_info = login_init_sid()
    # 获取任务列表
    ltm = ListTasksModel(num=max_num, username=login_info.username)
    aqm = AskQueueModel()
    layout = ListTaskLayout(ltm, aqm)
    layout.show()


class AskQueueModel:
    def __init__(self):
        self.num = None

    def ask(self):
        with UseTaskHttp() as http:
            queue_info = http.get("/task/queuing")
            self.num = queue_info["sum"]

    def table(self):
        qi = Table(
            expand=True,
            show_header=False,
            header_style="bold",
            title="[blue][b]Now Global Queue[/b]",
            highlight=True,
            border_style="blue",
        )
        qi.add_column("Queue Info", "Queue Info")
        self.ask()
        qi.add_row(f"[b]Task Queuing count: {self.num}[/b]")
        return qi


class ListTasksModel:
    def __init__(self, num: int, username: str):
        """
        :param num: 最大显示的任务数
        """
        self.num = num
        self.username = username

    def __dict__(self):
        return {"num": self.num}

    def list(self) -> List[TaskModel]:
        with UseTaskHttp() as http:
            tasks = http.get("/task", self.__dict__())
        return [TaskModel(self.username, task) for task in tasks]

    def table(self):
        st = Table(
            expand=True,
            show_header=True,
            title="[magenta][b]Now Task[/b]",
            highlight=True,
            border_style="magenta",
            show_lines=True,
        )
        st.add_column("Task ID", justify="right")
        st.add_column("Task Name", justify="center")
        st.add_column("Status", justify="center")
        st.add_column("URL", justify="center", no_wrap=True)
        st.add_column("Started Time", justify="center")
        st.add_column("Finished Time", justify="center")
        for tlm in self.list():
            status = tlm.status
            if status == "COMPLETED":
                status = f"[green]{status}[/green]"
            elif status == "CRASHED":
                status = f"[red]{status}[/red]"
            st.add_row(
                tlm.cuid,
                tlm.name,
                status,
                tlm.url,
                tlm.started_at,
                tlm.finished_at,
            )
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

    def __init__(self, ltm: ListTasksModel, aqm: AskQueueModel):
        self.layout = Layout()
        self.layout.split(Layout(name="header", size=3), Layout(name="main"))
        self.layout["main"].split_row(Layout(name="task_table", ratio=16), Layout(name="info_side", ratio=7))
        self.layout["info_side"].split_column(Layout(name="queue_info", ratio=1), Layout(name="datasets_list", ratio=5))
        self.layout["header"].update(ListTaskHeader())
        self.layout["task_table"].update(Panel(ltm.table(), border_style="magenta"))
        self.layout["queue_info"].update(Panel(aqm.table(), border_style="blue"))
        self.ltm = ltm
        self.aqm = aqm
        self.redraw_datasets_list()

    @property
    def datasets_list(self):
        to = Table(
            expand=True,
            show_header=True,
            header_style="bold",
            title="[blue][b]Datasets[/b]",
            highlight=True,
            show_lines=True,
            border_style="blue",
        )
        to.add_column("Dataset ID", justify="right")
        to.add_column("Dataset Name", justify="center")
        to.add_column("Dataset Desc", justify="center")
        to.add_column("Created Time", justify="center")
        return to

    def redraw_datasets_list(self):
        datasets_list = self.datasets_list
        with UseTaskHttp() as http:
            datasets = http.get("/task/datasets")
        for dataset in datasets:
            datasets_list.add_row(
                dataset["cuid"],
                dataset["name"],
                dataset.get("desc", ""),
                TaskModel.fmt_time(dataset["createdAt"]),
            )
        self.layout["datasets_list"].update(Panel(datasets_list, border_style="blue"))

    def show(self):
        with Live(self.layout, refresh_per_second=10, screen=True) as live:
            search_now = time.time()
            queue_now = time.time()
            while True:
                time.sleep(1)
                self.layout["header"].update(ListTaskHeader())
                if time.time() - search_now > 5:
                    search_now = time.time()
                    self.layout["task_table"].update(Panel(self.ltm.table(), border_style="magenta"))
                    self.redraw_datasets_list()
                if time.time() - queue_now > 3:
                    queue_now = time.time()
                    self.layout["queue_info"].update(Panel(self.aqm.table(), border_style="blue"))
                live.refresh()
