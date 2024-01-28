# 这里提供了 SwanLab 所有支持的图表类型

每个图表类型被抽象为一个组件,但是在外部通过本部分提供的[ChartsContainer](./ChartsContainer.vue)组件和[ChartContainer](./ChartsContainer.vue)组件统一调用,所需的就仅仅是传入chart配置和一些props

## 图表组件

对图表有以下要求:

每个图表组件有统一的 props，具体内容如下：

| props | 具体作用 |
| ----- | -------- |
| chart | 图表配置 |
| title | 图表标题 |

每个图表组件应该暴露一些通用 api，在其他地方通过 ref 获取组件对象进行调用，具体如下：

| api    | 具体作用                       |
| ------ | ------------------------------ |
| render | 渲染图表                       |
| change | 当图表内容改变时，重新渲染图表 |
| zoom   | 图表放大                       |

> 我们约定上述 api 有一个统一的传入参数 data，data 的数据结构为`{tag: tagData, ...}`  
> `tagData`是一个对象，每个key代表对应的数据，包含一个list，我们约定list的数据已经依据与元素内部的index完成了升序排序
> 具体的处理逻辑由各个图表组件自行决定和实现

大致的组件模版如下(js 部分), 可以输入`vue-chart`代码片段生成此模版(最新模版可能和下面的内容有所出入):

```js
import SLModal from '@swanlab-vue/components/SLModal.vue'
import * as UTILS from './utils'
import { ref,inject } from 'vue'

// ---------------------------------- 配置 ----------------------------------

const props = defineProps({
  title: {
    type: String,
    required: true
  },
  chart: {
    type: Object,
    required: true
  }
})

// 数据源 arrya
const source = props.chart.source
// 参考字段和显示名称
const { xField, xTitle } = UTILS.refrence2XField(props.chart.refrence)


// ---------------------------------- 图表颜色配置 ----------------------------------
// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const colors = inject('colors')
if (!colors) throw new Error('colors is not defined, please provide colors in parent component')

// ---------------------------------- 数据格式化 ----------------------------------
/**
 * 为了将数据格式化为图表可用的格式，需要将数据源中的数据进行格式化
 * @param { Object } data 待格式化的数据
 */
const format = (data) => {}

// ---------------------------------- 渲染功能 ----------------------------------
const render = (data) => {}

// ---------------------------------- 重渲染功能 ----------------------------------
const change = (data) => {}

// ---------------------------------- 放大功能 ----------------------------------
// 是否放大
const isZoom = ref(false)
// 放大数据
const zoom = (data) => {}

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  render,
  change,
  zoom
})
```
