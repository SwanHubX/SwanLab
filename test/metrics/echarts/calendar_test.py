"""
// @author: ComPleHN
// @file: calendar.vue
// @time: 2025/5/27 13:29
// @description: 本文件是对于echarts的 Calendar 图表测试 , 文件名不叫calendar是因为和库文件重名
"""
# ---------------------------------------------- Calendar - Calendar_heatmap ----------------------------------------------
import random
import datetime

import pyecharts.options as opts
from pyecharts.charts import Calendar


begin = datetime.date(2017, 1, 1)
end = datetime.date(2017, 12, 31)
data = [
    [str(begin + datetime.timedelta(days=i)), random.randint(1000, 25000)]
    for i in range((end - begin).days + 1)
]

c1 = (
    Calendar()
    .add(
        series_name="",
        yaxis_data=data,
        calendar_opts=opts.CalendarOpts(
            pos_top="120",
            pos_left="30",
            pos_right="30",
            range_="2017",
            yearlabel_opts=opts.CalendarYearLabelOpts(is_show=False),
        ),
    )
    .set_global_opts(
        title_opts=opts.TitleOpts(pos_top="30", pos_left="center", title="2017年步数情况"),
        visualmap_opts=opts.VisualMapOpts(
            max_=20000, min_=500, orient="horizontal", is_piecewise=False
        ),
    )
)

# ---------------------------------------------- Calendar - Calendar_label_setting ----------------------------------------------
import datetime
import random

from pyecharts import options as opts
from pyecharts.charts import Calendar


begin = datetime.date(2017, 1, 1)
end = datetime.date(2017, 12, 31)
data = [
    [str(begin + datetime.timedelta(days=i)), random.randint(1000, 25000)]
    for i in range((end - begin).days + 1)
]

c2 = (
    Calendar()
    .add(
        "",
        data,
        calendar_opts=opts.CalendarOpts(
            range_="2017",
            daylabel_opts=opts.CalendarDayLabelOpts(name_map="cn"),
            monthlabel_opts=opts.CalendarMonthLabelOpts(name_map="cn"),
        ),
    )
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Calendar-2017年微信步数情况(中文 Label)"),
        visualmap_opts=opts.VisualMapOpts(
            max_=20000,
            min_=500,
            orient="horizontal",
            is_piecewise=True,
            pos_top="230px",
            pos_left="100px",
        ),
    )
)

# ---------------------------------------------- Calendar - Calendar_base ----------------------------------------------
import datetime
import random

from pyecharts import options as opts
from pyecharts.charts import Calendar


begin = datetime.date(2017, 1, 1)
end = datetime.date(2017, 12, 31)
data = [
    [str(begin + datetime.timedelta(days=i)), random.randint(1000, 25000)]
    for i in range((end - begin).days + 1)
]

c3 = (
    Calendar()
    .add("", data, calendar_opts=opts.CalendarOpts(range_="2017"))
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Calendar-2017年微信步数情况"),
        visualmap_opts=opts.VisualMapOpts(
            max_=20000,
            min_=500,
            orient="horizontal",
            is_piecewise=True,
            pos_top="230px",
            pos_left="100px",
        ),
    )
)



import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="calendar",
    public=True,
)

swanlab.log(
    {
        "Calendar - Calendar_heatmap": c1,
        "Calendar - Calendar_label_setting": c2,
        "Calendar - Calendar_base": c3
    }
)



