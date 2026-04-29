import click
import orjson

from swanlab.api import Api
from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import format_output, save_output, with_custom_host
from swanlab.sdk.internal.pkg.safe import block


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
    default=None,
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
@with_custom_host
def create_user(username: str, password: str, api: Api):
    """Create a new user in the self-hosted instance."""
    err: list[str] = []
    resp: ApiResponseType | None = None
    with block(message=None, on_error=lambda e: err.append(str(e))):
        resp = api.self_hosted().create_user(username, password)
    format_output(resp if resp is not None else ApiResponseType(ok=False, errmsg=err[0]))


@selfhosted_cli.command("list-users")
@click.option("--page_num", "-n", default=1, type=int, help="Page number, default 1.")
@click.option("--page_size", "-s", default=20, type=int, help="Page size, default 20.")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all users.")
@with_custom_host
def list_users(page_num: int, page_size: int, fetch_all: bool, api: Api):
    """List users in the self-hosted instance."""
    err: list[str] = []
    users: list | None = None
    with block(message=None, on_error=lambda e: err.append(str(e))):
        users = list(api.self_hosted().get_users(page=page_num, size=page_size, all=fetch_all))
    if users is not None:
        format_output(ApiResponseType(ok=True, data={"list": users}))
    else:
        format_output(ApiResponseType(ok=False, errmsg=err[0]))
