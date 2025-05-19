import swanlab

bar = (
    swanlab.echarts.Bar()
    .add_xaxis(["衬衫", "羊毛衫", "雪纺衫"])
    .add_yaxis("商家A", [5, 20, 36])
)

swanlab.log({"echart.bar": bar})