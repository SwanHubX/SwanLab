"""
@author: cunyue
@file: video.py
@time: 2025/7/17 11:32
@description: 视频 gif 测试 demo
"""

import os.path
import random

from PIL import Image as PILImage
from PIL import ImageDraw

import swanlab
from tutils import TEMP_PATH

swanlab.init()


gif_path = os.path.join(TEMP_PATH, "test.gif")


# 创建一个GIF动画
def create_mock_gif(output_path, width=200, height=200, frames=10, duration=100):
    """
    创建一个简单的mock GIF动画

    参数:
        output_path: 输出GIF文件路径
        width: 图像宽度(像素)
        height: 图像高度(像素)
        frames: 动画帧数
        duration: 每帧显示时间(毫秒)
    """
    images = []

    for i in range(frames):
        # 创建一个新的RGB图像
        img = PILImage.new('RGB', (width, height), color=(255, 255, 255))
        draw = ImageDraw.Draw(img)

        # 随机生成颜色
        r = random.randint(0, 255)
        g = random.randint(0, 255)
        b = random.randint(0, 255)

        # 在图像上绘制一个随机的圆形
        x = random.randint(0, width)
        y = random.randint(0, height)
        radius = random.randint(10, min(width, height) // 2)
        draw.ellipse([x - radius, y - radius, x + radius, y + radius], fill=(r, g, b))

        # 将当前帧添加到列表中
        images.append(img)

    # 保存为GIF动画
    images[0].save(output_path, save_all=True, append_images=images[1:], duration=duration, loop=0)


create_mock_gif(gif_path, width=300, height=300, frames=15, duration=200)

swanlab.log({"video": swanlab.Video(gif_path, caption="This is a test video")})
