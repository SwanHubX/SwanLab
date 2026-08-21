from __future__ import annotations

from typing import Optional

import click

from swanlab.api import Api
from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import api_command


@click.group("self-hosted")
def selfhosted_cli():
    """Self-Hosted deployment management commands."""
    pass


@selfhosted_cli.command("info")
@api_command
def get_info(api: Api) -> ApiResponseType:
    """Show self-hosted instance info."""
    return api.self_hosted().wrapper()


@selfhosted_cli.command("create-user")
@click.option("--username", "-u", type=str, required=True, help="username to create")
@click.option("--password", "-p", type=str, required=True, help="password to create")
@api_command
def create_user(username: str, password: str, api: Api) -> ApiResponseType:
    """Create a new user in the self-hosted instance."""
    return api.self_hosted().create_user(username, password)


@selfhosted_cli.command("list-users")
@click.option("--page_num", "-n", default=1, type=int, help="Page number, default 1.")
@click.option("--page_size", "-s", default=20, type=int, help="Page size, default 20.")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all users.")
@api_command
def list_users(page_num: int, page_size: int, fetch_all: bool, api: Api) -> ApiResponseType:
    """List users in the self-hosted instance."""
    sh = api.self_hosted()
    users = list(sh.get_users(page=page_num, size=page_size, all=fetch_all))
    if sh._errors:
        return ApiResponseType(ok=False, errmsg="; ".join(sh._errors))
    return ApiResponseType(ok=True, data={"list": users})


@selfhosted_cli.command("list-projects")
@click.option("--page_num", "-n", default=1, type=int, help="Page number, default 1.")
@click.option("--page_size", "-s", default=20, type=int, help="Page size, default 20.")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all projects.")
@click.option("--search", default=None, type=str, help="Search keyword.")
@click.option("--creator", default=None, type=str, help="Filter by creator username.")
@click.option("--workspace", default=None, type=str, help="Filter by workspace username.")
@api_command
def list_projects(
    page_num: int,
    page_size: int,
    fetch_all: bool,
    search: Optional[str],
    creator: Optional[str],
    workspace: Optional[str],
    api: Api,
) -> ApiResponseType:
    """List all projects in the self-hosted instance."""
    sh = api.self_hosted()
    projects = list(
        sh.get_projects(
            page=page_num,
            size=page_size,
            all=fetch_all,
            search=search,
            creator=creator,
            group=workspace,
        )
    )
    if sh._errors:
        return ApiResponseType(ok=False, errmsg="; ".join(sh._errors))
    return ApiResponseType(ok=True, data={"list": projects})


@selfhosted_cli.command("summary")
@api_command
def get_summary(api: Api) -> ApiResponseType:
    """Show system usage summary (root only)."""
    return api.self_hosted().get_usage_summary()


@selfhosted_cli.command("list-workspaces")
@click.option("--page_num", "-n", default=1, type=int, help="Page number, default 1.")
@click.option("--page_size", "-s", default=20, type=int, help="Page size, default 20.")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all workspaces.")
@click.option("--search", default=None, type=str, help="Search keyword.")
@api_command
def list_groups(page_num: int, page_size: int, fetch_all: bool, search: Optional[str], api: Api) -> ApiResponseType:
    """List all workspaces in the self-hosted instance."""
    sh = api.self_hosted()
    workspaces = list(
        sh.get_groups(
            page=page_num,
            size=page_size,
            all=fetch_all,
            search=search,
        )
    )
    if sh._errors:
        return ApiResponseType(ok=False, errmsg="; ".join(sh._errors))
    return ApiResponseType(ok=True, data={"list": workspaces})
