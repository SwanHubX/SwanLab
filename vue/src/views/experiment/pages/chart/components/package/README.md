# 这里提供了 SwanLab 所有支持的图表类型

每个图表类型被抽象为一个组件，在 [ChartContainer](../ChartContainer.vue) 中选择性渲染。

## 图表组件

对图表有以下要求:

每个图表组件有统一的 props，具体内容如下：

| props | 具体作用                |
| ----- | ----------------------- |
| chart | 图表配置                |
| data  | 基于 sources 的图表数据 |

每个图表组件应该暴露一些通用 api，在[ChartContainer](../ChartContainer.vue)中通过 ref 获取组件对象进行调用，具体如下：

| api    | 具体作用                       |
| ------ | ------------------------------ |
| render | 渲染图表                       |
| change | 当图表内容改变时，重新渲染图表 |
| zoom   | 图表放大                       |
