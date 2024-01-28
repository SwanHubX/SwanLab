# 这里提供了 SwanLab 所有支持的图表类型

每个图表类型被抽象为一个组件,但是在外部通过本部分提供的[ChartsContainer](./ChartsContainer.vue)组件统一调用，所需的就是传入chart配置  
一个ChartsContainer组件代表一个命名空间

## 图表组件

所有的图表组件都将有统一的props输入：
| 属性 | 类型 | 说明 |
| --- | --- | --- |
| title | String | 图表标题 |
| chart | Object | 图表配置 |
| index | Number | 图表在当前命名空间中的排序 |

所有的图表也将有通用的api，通过defineExpose暴露给父组件：
| 方法 | 类型 | 说明 |
| --- | --- | --- |
| render | Function | 渲染图表 |
| change | Function | 重渲染图表 |
| zoom | Function | 放大图表 |

其中，这些通用的api都将接受一个参数，即图表数据，其数据格式为：

```js
const data = {
    "数据名称": {
        list:[
            data: "数据",
            step: "步骤"，
            create_time: "创建时间"
        ],
        ... // 一些其他的配置字段
    }
    ...
}
```

本模块不再遵守vue单向数据流规范，但是也不鼓励暴露数据，而应该暴露api，通过api来完成逻辑的处理

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
