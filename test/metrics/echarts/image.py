"""
// @author: ComPleHN
// @file: image.vue
// @time: 2025/5/27 17:04
// @description: 本文件是对于echarts的 图片 的测试
"""
# ---------------------------------------------- Image - Image_base ----------------------------------------------
from pyecharts.components import Image
from pyecharts.options import ComponentTitleOpts


image = Image()

img_src = (
    "https://user-images.githubusercontent.com/19553554/"
    "71825144-2d568180-30d6-11ea-8ee0-63c849cfd934.png"
)
image.add(
    src=img_src,
    style_opts={"width": "200px", "height": "200px", "style": "margin-top: 20px"},
)
image.set_global_opts(
    title_opts=ComponentTitleOpts(title="Image-基本示例", subtitle="我是副标题支持换行哦")
)

c1 = image

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="image",
    public=True,
)

swanlab.log(
    {
        "Calendar - Calendar_heatmap": c1,
    }
)
