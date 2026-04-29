from __future__ import annotations

from typing import Optional

import click
import orjson

from swanlab.api import Api
from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import format_output, save_output, with_custom_host


@click.group("selfhosted")
def selfhosted_cli():
    """Self-hosted deployment management commands."""
    pass


@selfhosted_cli.command("info")
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    help="Save output as JSON to current directory.",
)
@with_custom_host
def get_info(save_name: str, api: Api):
    """Show self-hosted instance info."""
    resp = api.self_hosted().wrapper()
    payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@selfhosted_cli.command("create-user")
@click.option("--username", "-u", type=str, required=True, help="username to create")
@click.option("--password", "-p", type=str, required=True, help="password to create")
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    help="Save output as JSON to current directory.",
)
@with_custom_host
def create_user(username: str, password: str, save_name: str, api: Api):
    """Create a new user in the self-hosted instance."""
    try:
        resp = api.self_hosted().create_user(username, password)
    except ValueError as e:
        resp = ApiResponseType(ok=False, errmsg=str(e))
    payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@selfhosted_cli.command("list-users")
@click.option("--page_num", "-n", default=1, type=int, help="Page number, default 1.")
@click.option("--page_size", "-s", default=20, type=int, help="Page size, default 20.")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all users.")
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    help="Save output as JSON to current directory.",
)
@with_custom_host
def list_users(page_num: int, page_size: int, fetch_all: bool, save_name: str, api: Api):
    """List users in the self-hosted instance."""
    try:
        users = list(api.self_hosted().get_users(page=page_num, size=page_size, all=fetch_all))
    except ValueError as e:
        payload = format_output(ApiResponseType(ok=False, errmsg=str(e)))
    else:
        resp = ApiResponseType(ok=True, data={"list": users})
        payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@selfhosted_cli.command("list-projects")
@click.option("--page_num", "-n", default=1, type=int, help="Page number, default 1.")
@click.option("--page_size", "-s", default=20, type=int, help="Page size, default 20.")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all projects.")
@click.option("--search", default=None, type=str, help="Search keyword.")
@click.option("--creator", default=None, type=str, help="Filter by creator username.")
@click.option("--workspace", default=None, type=str, help="Filter by workspace username.")
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    help="Save output as JSON to current directory.",
)
@with_custom_host
def list_projects(
    page_num: int,
    page_size: int,
    fetch_all: bool,
    search: Optional[str],
    creator: Optional[str],
    workspace: Optional[str],
    save_name: str,
    api: Api,
):
    """List all projects in the self-hosted instance."""
    try:
        projects = list(
            api.self_hosted().get_projects(
                page=page_num,
                size=page_size,
                all=fetch_all,
                search=search,
                creator=creator,
                group=workspace,
            )
        )
    except ValueError as e:
        payload = format_output(ApiResponseType(ok=False, errmsg=str(e)))
    else:
        resp = ApiResponseType(ok=True, data={"list": projects})
        payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@selfhosted_cli.command("summary")
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    help="Save output as JSON to current directory.",
)
@with_custom_host
def get_summary(save_name: str, api: Api):
    """Show system usage summary (root only)."""
    try:
        resp = api.self_hosted().get_usage_summary()
    except ValueError as e:
        resp = ApiResponseType(ok=False, errmsg=str(e))
    payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@selfhosted_cli.command("list-workspaces")
@click.option("--page_num", "-n", default=1, type=int, help="Page number, default 1.")
@click.option("--page_size", "-s", default=20, type=int, help="Page size, default 20.")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all workspaces.")
@click.option("--search", default=None, type=str, help="Search keyword.")
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    help="Save output as JSON to current directory.",
)
@with_custom_host
def list_groups(page_num: int, page_size: int, fetch_all: bool, search: Optional[str], save_name: str, api: Api):
    """List all workspaces in the self-hosted instance."""
    try:
        workspaces = list(
            api.self_hosted().get_groups(
                page=page_num,
                size=page_size,
                all=fetch_all,
                search=search,
            )
        )
    except ValueError as e:
        payload = format_output(ApiResponseType(ok=False, errmsg=str(e)))
    else:
        resp = ApiResponseType(ok=True, data={"list": workspaces})
        payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)
