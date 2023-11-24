from fastapi import FastAPI
from fastapi.responses import HTMLResponse, FileResponse
from fastapi.staticfiles import StaticFiles
import uvicorn
import threading
import os

"""
在此处注册静态文件路径，因为静态文件由vue框架编译后生成，在配置中，编译后的文件存储在/swanlab/template中
入口文件为index.html，网页图标为logo.ico，其他文件为assets文件夹中的文件
因此，需要指定文件路径与文件名，用于后端服务的响应（这在下面的路由配置中也有说明）
"""
FILEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE_PATH = os.path.join(FILEPATH, "template")
INDEX = os.path.join(TEMPLATE_PATH, "index.html")
ASSETS = os.path.join(TEMPLATE_PATH, "assets")


# 服务全局变量
_app = FastAPI()


# 响应首页
@_app.get("/", response_class=HTMLResponse)
async def _():
    # 读取 HTML 文件内容并返回
    with open(INDEX, "r", encoding="utf-8") as file:
        html_content = file.read()
    return HTMLResponse(content=html_content, status_code=200)


# 响应logo内容
@_app.get("/swanlab.svg")
async def _():
    return FileResponse(os.path.join(TEMPLATE_PATH, "swanlab.svg"))


# 响应其他assets文件，FileResponse会自动根据文件类型设置响应头
@_app.get("/assets/{file_name}")
async def _(file_name: str):
    path = os.path.join(ASSETS, file_name)
    return FileResponse(path)


class SwanWeb(object):
    """SwanWeb类用于启动函数，用于连接前端与interface层，实现服务的开启和数据的传输
    此类将被SwanLab继承
    """

    # 是否正在运行
    def __init__(self, config={}):
        """
        初始化服务，原则上这里的参数不应该被使用本库的人所修改，因为这里是配置一些底层参数的，与用户无关
        有额外的init方法给用户配置相关所需参数
        """
        self.logs = {}  # 初始化为一个字典，存放日志
        self.server_thread = None  # 服务线程
        self.config = config  # 配置参数

    def run(self):
        """启动服务，执行此函数，将web服务开启在子线程
        开启web服务
        """
        static = StaticFiles(directory=ASSETS)
        _app.mount(ASSETS, static, name="assert")
        # 在初始化时就启动服务器
        self.server_thread = threading.Thread(target=self.__run_server)
        self.server_thread.start()

    # 开启服务
    def __run_server(self):
        uvicorn.run(_app, host="127.0.0.1", port=8080, log_level="info")

    # 生成log
    def log(self, log):
        self.logs.setdefault(log["name"], []).append(log["data"])
