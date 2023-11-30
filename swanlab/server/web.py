from fastapi.staticfiles import StaticFiles
import uvicorn
import threading
import mimetypes
import os
from ..database import SwanDataBase
from ..utils import color
from .route import app as _app

"""
在此处注册静态文件路径，因为静态文件由vue框架编译后生成，在配置中，编译后的文件存储在/swanlab/template中
入口文件为index.html，网页图标为logo.ico，其他文件为assets文件夹中的文件
因此，需要指定文件路径与文件名，用于后端服务的响应（这在下面的路由配置中也有说明）
"""
# 注册静态文件路径
mimetypes.add_type("application/javascript", ".js")
mimetypes.add_type("text/css", ".css")
FILEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE_PATH = os.path.join(FILEPATH, "template")
ASSETS = os.path.join(TEMPLATE_PATH, "assets")


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
        print("SwanLab server is running on " + color.green("http://{}:{}").format(host, self.port))
        if not self.share:
            print("You can share this server by setting share=True")

    def init(self, log_path: str, port: int = 10101, share: bool = False, info_level: str = "info"):
        """初始化服务，用于配置服务参数以及启动服务，加载数据库

        Parameters
        ----------
        log_path : str
            日志文件夹路径，运行的日志将会保存在此文件夹中
        port : int, optional
            服务开启的端口, by default 10101
        share : bool, optional
            是否需要开启局域网共享，这意味着将跑在localhost上, by default False
        info_level : str, optional
            日志等级，可选参数为: debug,info,warning,error,critical, by default "info"
        """
        self.share = share
        self.port = port
        self.info_level = info_level
        self.__run()
