from fastapi import FastAPI
from fastapi.responses import HTMLResponse, FileResponse, JSONResponse
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


# 开启服务
def run_server(host, port, log_level="warning"):
    """开启服务
    此函数将在子线程中运行
    """
    log_level = "info"
    uvicorn.run(_app, host=host, port=port, log_level=log_level)


class SwanWeb(object):
    """SwanWeb类用于启动函数，用于连接前端与interface层，实现服务的开启和数据的传输
    此类将被SwanLab继承
    """

    # 是否正在运行
    def __init__(self, share=False):
        """
        初始化服务，原则上这里的参数不应该被使用本库的人所修改，因为这里是配置一些底层参数的，与用户无关
        有额外的init方法给用户配置相关所需参数
        """
        self.logs = {}  # 初始化为一个字典，存放日志
        self.server_thread = None  # 服务线程
        self.share = share  # 是否开启当前服务网络共享，默认为False，代表当前服务跑在127.0.0.1上，只能本机访问
        self.port = 10101

    def run(self):
        """启动服务，执行此函数，将web服务开启
        由于用到了子线程，因此服务实际上启动在子线程中，主线程将继续执行
        """
        static = StaticFiles(directory=ASSETS)
        # 将assets文件夹注册为静态文件路径，这样不再需要单独响应每个文件
        _app.mount("/assets", static, name="assets")
        host = "127.0.0.1" if not self.share else "0.0.0.0"
        self.server_thread = threading.Thread(target=run_server, args=(host, self.port))
        self.server_thread.start()
        # 日志打印
        print("SwanLab server is running on http://{}:{}".format(host, self.port))
        if not self.share:
            print("You can share this server by setting share=True")

    # 生成log
    def log(self, log):
        self.logs.setdefault(log["name"], []).append(log["data"])
