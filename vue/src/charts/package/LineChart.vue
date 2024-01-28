<template>
  <!-- 图表标题 -->
  <p class="text-center font-semibold">{{ title }}</p>
  <div class="flex flex-col justify-center grow text-dimmer gap-2" v-if="error">
    <SLIcon class="mx-auto h-5 w-5" icon="error" />
    <p class="text-center text-sm">
      {{ $t('experiment.chart.charts.line.error', { type: error['data_class'], tag: source[0] }) }}
    </p>
  </div>
  <template v-else>
    <!-- x轴坐标单位 -->
    <p class="absolute right-5 bottom-10 text-xs text-dimmer scale-90">{{ xTitle }}</p>
    <!-- 图表主体 -->
    <div ref="g2Ref"></div>
    <!-- 放大效果 -->
    <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom">
      <p class="text-center mt-4 mb-10 text-2xl font-semibold">{{ title }}</p>
      <div ref="g2ZoomRef"></div>
      <p class="absolute right-12 bottom-16 text-xs text-dimmer scale-90">{{ xTitle }}</p>
    </SLModal>
  </template>
</template>

<script setup>
/**
 * @description: 折线图表
 * @file: LineChart.vue
 * @since: 2023-12-25 20:17:19
 **/
import SLModal from '@swanlab-vue/components/SLModal.vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import { Line } from '@antv/g2plot'
import * as UTILS from './utils'
import { ref, inject } from 'vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import { formatNumber2SN } from '@swanlab-vue/utils/common'
import { useProjectStore } from '@swanlab-vue/store'

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

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------

const error = ref(props.chart.error)

// ---------------------------------- 图表颜色配置 ----------------------------------
// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const colors = inject('colors')
if (!colors) throw new Error('colors is not defined, please provide colors in parent component')

// ---------------------------------- 组件渲染逻辑 ----------------------------------

// 组件对象
const g2Ref = ref()
const g2ZoomRef = ref()
// 数据源 arrya
const source = props.chart.source
// 参考字段和显示名称
const reference = props.chart.reference
// 拿到参考系，未来图表可能有不同的x轴依据，比如step、time等，这里需要根据设置的reference来决定
const { xField, xTitle } = UTILS.refrence2XField[reference]
// 默认y轴的依据key是data
const yField = 'data'
const seriesField = 'series'
const colorField = 'type'
// 创建图表的函数
const createChart = (dom, data, config = { interactions: undefined, height: 200 }) => {
  const c = new Line(dom, {
    data,
    // 默认的x轴依据key为step
    xField,
    // 默认的y轴依据key为data
    yField,
    // 多数据的时候，需要设置seriesField，单数据也可以设置，但是没什么区别
    seriesField,
    colorField,
    color: colors,
    // 坐标轴相关
    xAxis: {
      // 自定义坐标轴的刻度，暂时没有找到文档，通过源码来看是返回一个数组，数组内是字符串，代表刻度
      tickMethod: formatXAxisTick,
      // 在此处完成X轴数据的格式化
      label: {
        formatter: (data) => {
          return formatNumber2K(data)
        }
      }
    },
    yAxis: {
      min: null,
      label: {
        // 在此处完成Y轴数据的格式化
        formatter: (data) => {
          return formatNumber2SN(data)
        }
      },
      // 显示y轴刻度线#bfbfbf
      line: {
        style: {
          stroke: '#bfbfbf'
        }
      },
      tickLine: {
        style: {
          stroke: '#bfbfbf'
        }
      }
    },
    tooltip: {
      // 在此处完成悬浮数据提示的格式化
      // 如果需要自定义浮窗，可以用下面的customContent，但是目前不管
      formatter: (data) => {
        // console.log(data)
        return { name: data.series, value: formatNumber2SN(data.data) }
      }
      // customContent: (title, data) => {
      //   console.log(title, data)
      //   return `<div>${title}</div>`
      // }
    },
    // 大小相关
    height: 200,
    width: undefined,
    autoFit: true,
    // 开启一些交互
    interactions: undefined,
    // 样式相关
    // smooth: true, // 平滑曲线
    // 线条样式
    // lineStyle,
    ...config
  })
  c.render()
  return c
}

// ---------------------------------- 数据格式化 ----------------------------------
/**
 * 为了将数据格式化为图表可用的格式，需要将数据源中的数据进行格式化
 * 遍历data的所有key，合并其中的list为一个数组
 * @param { Object } data 待格式化的数据
 * @returns { Object } 格式化后的数据, { d: [{}, {}, ...], config: {} } config是图表的一些其他配置
 */
const format = (data) => {
  // 如果source的长度小于1，抛出错误
  if (source.length < 1) throw new Error('source length must be greater than 1')
  // 新的数据,遍历得到
  const d = []
  Object.keys(data).forEach((key) => {
    // 如果不是单数据，需要将所有数据的list合并为一个数组
    data[key].list.forEach((item) => {
      // item新加series字段，用于标识数据来源
      d.push({ ...item, series: key })
    })
  })
  return { d }
}

/**
 * 以千为单位格式化数字，例如:
 * 100 => 100 (如果不是1000的倍数，则直接返回)
 * 1000 => 1k
 * 10000 => 10k
 * @param {number} num 待格式化的数字
 * @returns {string} 格式化后的字符串
 */
const formatNumber2K = (num) => {
  if (num % 1000 !== 0 || num == 0) return String(num)
  return `${num / 1000}k`
}

/**
 * 格式化x轴的刻度，最终返回一个数组，数组内是字符串，代表刻度
 * 所有的处理都基于category的values，即category.values是一个数组，数组内是字符串，代表所有数据的刻度，为字符串类型，已经基于数值大小完成了升序排序
 * 划分规则如下：
 * 1. [0, 8)，直接返回
 * 2. [8, 15)，每隔2个取一个（第一个和最后一个都取） 15=2*7+1
 * 3. [15, 35), 每隔5个取一个（第一个和最后一个都取） 35=5*7
 * 4. [35, 70), 每隔10个取一个（第一个和最后一个都取） 70=10*7
 * 5. [70, 140), 每隔20个取一个（第一个和最后一个都取） 140=20*7
 * 6. [140, 350), 每隔50个取一个（第一个和最后一个都取） 350=50*7
 * ...
 *
 * @param {object} category 目录配置
 * @returns {array} 刻度数组
 */
const formatXAxisTick = (category) => {
  // console.log('category', category)
  const { values } = category
  const length = values.length
  if (length < 8) return values
  if (length < 15) return values.filter((_, i) => i % 2 === 0)
  if (length <= 35) return values.filter((_, i) => i % 5 === 0)
  // 接下来，将length归一化到[35, 350)区间内，判断步长取值
  // 1. 计算归一化尺度，这应该是10的n次方,具体计算是/35，然后判断其小于哪个10的n次方
  const p = Math.ceil(Math.log10(length / 35)) - 1
  // 现在区间就确定了，是[35*10**(p-1), 35*10**p]
  // 然后判断一下在哪个区间内，这样就可以取步长
  const lengthScale = length / 10 ** p
  let step
  if (lengthScale < 70) {
    step = 10 * 10 ** p
  } else if (lengthScale < 140) {
    step = 20 * 10 ** p
  } else {
    step = 50 * 10 ** p
  }
  // console.log('p', p)
  // console.log('step', step)
  // 判断这个长度能容纳多少个步长
  const count = Math.floor(length / step)
  // console.log('count', count)
  // 生成刻度，与步长绑定,去重
  const ticks = new Set()
  // 以最小值为基准，生成刻度
  const base = Math.floor(Number(values[0]) / step)
  // 刻度不大于最大值
  const max = Math.floor(Number(values[length - 1]))
  // console.log('base', base)
  for (let i = base; i <= base + count; i++) {
    i * step < max && ticks.add(i * step)
  }
  // set转array
  console.log('ticks', Array.from(ticks))
  return Array.from(ticks)
}

// ---------------------------------- 渲染、重渲染功能 ----------------------------------
let chartObj = null
// 渲染
const render = (data) => {
  // console.log('渲染折线图')
  // console.log('data', data)
  const { d } = format(data)
  // console.log('data', data)
  chartObj = createChart(g2Ref.value, d)
  // console.log('chartObj', chartObj)
  // 可以使用update api来更新配置
}
// 重渲染
const change = (data) => {
  const { d } = format(data)
  // 如果d的长度超过1200, change函数等于render函数
  if (d.length > 1200) {
    chartObj.destroy()
    return render(data)
  }

  // console.log('更新...')
  // console.log(chartObj)
  // updateYAxis(yAxis)
  chartObj.changeData(d)
}

// ---------------------------------- 放大功能 ----------------------------------
// 是否放大
const isZoom = ref(false)
// 放大数据
const zoom = (data) => {
  isZoom.value = true
  // 当前window的高度
  const { d } = format(data)
  const height = window.innerHeight * 0.6
  addTaskToBrowserMainThread(() => {
    createChart(g2ZoomRef.value, d, {
      interactions: [{ type: 'brush-x' }],
      height
    })
  })
}

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  render,
  change,
  zoom
})
</script>

<style lang="scss" scoped></style>
