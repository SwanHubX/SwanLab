from fastapi import FastAPI
from fastapi.responses import HTMLResponse, FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from typing import Union
import uvicorn
import threading
import os
from ..database import SwanDataBase

"""
在此处注册静态文件路径，因为静态文件由vue框架编译后生成，在配置中，编译后的文件存储在/swanlab/template中
入口文件为index.html，网页图标为logo.ico，其他文件为assets文件夹中的文件
因此，需要指定文件路径与文件名，用于后端服务的响应（这在下面的路由配置中也有说明）
"""
FILEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE_PATH = os.path.join(FILEPATH, "template")
INDEX = os.path.join(TEMPLATE_PATH, "index.html")
ASSETS = os.path.join(TEMPLATE_PATH, "assets")

# 服务全局对象
_app = FastAPI()


# 响应首页
@_app.get("/", response_class=HTMLResponse)
async def _():
    # 读取 HTML 文件内容并返回
    with open(INDEX, "r", encoding="utf-8") as file:
        html_content = file.read()
    return HTMLResponse(content=html_content, status_code=200)


# 响应logo内容
# TODO 后续可以考虑将logo.ico放在assets中，这样就不需要单独响应了
@_app.get("/logo.ico")
async def _():
    return FileResponse(os.path.join(TEMPLATE_PATH, "logo.ico"))


import random


# 测试路由，每次请求返回一个0到30的随机数
@_app.get("/api/test")
async def _():
    # 生成一个 0 到 30 之间的随机整数
    random_number = random.randint(0, 30)
    return JSONResponse({"data": random_number}, status_code=200, headers={"Access-Control-Allow-Origin": "*"})


class SwanWeb(object):
    """SwanWeb类用于启动函数，用于连接前端与interface层，实现服务的开启和数据的传输
    此类将被SwanLab继承
    """

    # 是否正在运行
    def __init__(self):
        """
        初始化服务，原则上这里的参数不应该被使用本库的人所修改，因为这里是配置一些底层参数的，与用户无关
        有额外的init方法给用户配置相关所需参数
        """
        self.server_thread = None  # 服务线程
        self.database: SwanDataBase = SwanDataBase()  # 数据库对象

        # 可配置内容
        self.share = False  # 是否开启当前服务网络共享，默认为False，代表当前服务跑在127.0.0.1上，只能本机访问
        self.port = 10101
        self.info_level = "info"

    def __run(self):
        """启动服务，执行此函数，将web服务开启
        由于用到了子线程，因此服务实际上启动在子线程中，主线程将继续执行
        """

        # 开启服务的函数
        def run_server(host, port, log_level="warning"):
            """开启服务
            此函数将在子线程中运行
            """
            uvicorn.run(_app, host=host, port=port, log_level=log_level)

        # 注册静态文件
        static = StaticFiles(directory=ASSETS)
        # 将assets文件夹注册为静态文件路径，这样不再需要单独响应每个文件
        _app.mount("/assets", static, name="assets")
        host = "127.0.0.1" if not self.share else "localhost"
        self.server_thread = threading.Thread(target=run_server, args=(host, self.port))
        self.server_thread.start()
        # 日志打印
        print("SwanLab server is running on http://{}:{}".format(host, self.port))
        if not self.share:
            print("You can share this server by setting share=True")

    def init(self, port: int = 10101, share: bool = False, info_level: str = "info"):
        """初始化服务，用于配置服务参数以及启动服务

        Parameters
        ----------
        port : int, optional
            服务开启的端口, by default 10101
        share : bool, optional
            是否需要开启局域网共享，这意味着将跑在localhost上, by default False
        info_level : str, optional
            日志等级，可选参数为：debug,info,warning,error,critical, by default "info"
        """
        self.share = share
        self.port = port
        self.info_level = info_level
        self.__run()

    # 生成log
    def log(self, tag: str, index: int, data: Union[int, float], namspace: str = "default"):
        """生成日志，用于在前端显示

        Parameters
        ----------
        tag : str
            添加数据的标签，用于区分不同的数据
        index : int
            数据的索引，用于区分不同的数据
        data : Union[int, float]
            添加的数据，暂时只支持int和float类型
        namspace : str, optional
            命名空间，用于区分不同的数据，system用于存储系统信息， by default "default"
        """
        self.database.add(tag, index, data, namspace)
