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
from .utils import TaskModel, UseTaskHttp, login_init_sid


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
    # 获取任务列表
    ltm = ListTasksModel(num=max_num, username=login_info.username)
    layout = ListTaskLayout(ltm)
    layout.show()


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

    def __init__(self, ltm: ListTasksModel):
        self.event = []
        self.add_event(f"👏Welcome, [b]{ltm.username}[/b].")
        self.add_event("⌛️Task board is loading...")
        self.layout = Layout()
        self.layout.split(
            Layout(name="header", size=3),
            Layout(name="main")
        )
        self.layout["main"].split_row(
            Layout(name="task_table", ratio=4),
            Layout(name="term_output", ratio=1)
        )
        self.layout["header"].update(ListTaskHeader())
        self.layout["task_table"].update(Panel(ltm.table(), border_style="magenta"))
        self.ltm = ltm
        self.add_event("🍺Task board is loaded.")
        self.redraw_term_output()

    @property
    def term_output(self):
        to = Table(
            expand=True,
            show_header=False,
            header_style="bold",
            title="[blue][b]Log Messages[/b]",
            highlight=True,
            border_style="blue",
            show_footer=True,
            footer_style="bold",
        )
        to.add_column(
            "Log Output",
            "Run [b][white]swanlab task search [Task ID][/white][/b] to get more task info"
        )
        return to

    def redraw_term_output(self, ):
        term_output = self.term_output
        for row in self.event:
            term_output.add_row(row)
        self.layout["term_output"].update(Panel(term_output, border_style="blue"))

    def add_event(self, info: str, max_length=15):
        # 事件格式：yyyy-mm-dd hh:mm:ss - info
        self.event.append(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {info}")
        while len(self.event) > max_length:
            self.event.pop(0)

    def show(self):
        with Live(self.layout, refresh_per_second=10, screen=True) as live:
            now = time.time()
            while True:
                time.sleep(1)
                self.layout["header"].update(ListTaskHeader())
                if time.time() - now > 5:
                    now = time.time()
                    self.add_event("🔍Searching for new tasks...")
                    self.layout["task_table"].update(Panel(self.ltm.table(), border_style="magenta"))
                    self.redraw_term_output()
                live.refresh()
