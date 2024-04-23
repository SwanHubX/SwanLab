from IPython.core.magic import Magics, magics_class, line_magic, line_cell_magic
from IPython.core.magic_arguments import argument, magic_arguments, parse_argstring
import subprocess

#  嵌入SwanLab网页
# class IFrame:
#     def __init__(self, path=None, opts=None):
#         self.path = path
#         self.opts = opts or {}
#         self.displayed = False
#         self.height = self.opts.get("height", 420)

#     def maybe_display(self) -> bool:
#         if not self.displayed and (self.path or wandb.run):
#             display(self)
#         return self.displayed

#     def _repr_html_(self):
#         try:
#             self.displayed = True
#             if self.opts.get("workspace", False):
#                 if self.path is None and wandb.run:
#                     self.path = wandb.run.path
#             if isinstance(self.path, str):
#                 object = self.api.from_path(self.path)
#             else:
#                 object = wandb.run
#             if object is None:
#                 if wandb.Api().api_key is None:
#                     return "You must be logged in to render wandb in jupyter, run `wandb.login()`"
#                 else:
#                     object = self.api.project(
#                         "/".join(
#                             [
#                                 wandb.Api().default_entity,
#                                 wandb.util.auto_project_name(None),
#                             ]
#                         )
#                     )
#             return object.to_html(self.height, hidden=False)
#         except wandb.Error as e:
#             return f"Can't display wandb interface<br/>{e}"


@magics_class
class SwanLabMagics(Magics):
    def __init__(self, shell, require_interaction=False):
        super().__init__(shell)
        self.options = {}

    @magic_arguments()
    @argument(
        "-l",
        "--logdir",
        default="swanlog",
        help="指定swanlab的日志文件位置",
    )
    @line_magic
    def swanlab(self, line, cell=None):
        args = parse_argstring(self.swanlab, line)
        logdir = args.logdir

        result = subprocess.run(
            [
                "swanlab",
                "watch",
                "-l",
                f"{logdir}",
            ],
            shell=True,
            text=True,
            capture_output=True,
        )

        if result.returncode != 0:
            print("Error: Command failed")
            print(result.stderr)
            return result.stderr

        return result.stdout

        # if line:
        #     return f"Hello logdir {logdir}!"
