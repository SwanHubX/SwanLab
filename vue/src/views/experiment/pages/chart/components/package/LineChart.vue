<template>
  <!-- 图表标题 -->
  <p class="text-center font-semibold">{{ title }}</p>
  <div class="flex flex-col justify-center grow text-dimmer gap-2" v-if="error">
    <SLIcon class="mx-auto h-5 w-5" icon="error" />
    <p class="text-center text-xs">
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
import { ref } from 'vue'
import { useExperimentStroe } from '@swanlab-vue/store'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'

// ---------------------------------- 配置 ----------------------------------
const experimentStore = useExperimentStroe()
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

// ---------------------------------- 组件渲染逻辑 ----------------------------------

// 组件对象
const g2Ref = ref()
const g2ZoomRef = ref()
// 数据源 arrya
const source = props.chart.source
// 参考字段和显示名称
const reference = props.chart.reference
// 图表默认颜色
const defaultColor = props.chart.config?.color || experimentStore.defaultColor
// 拿到参考系
const { xField, xTitle } = UTILS.refrence2XField[reference]
// 创建图表的函数
// FIXME 兼容多数据情况
const createChart = (dom, data, config = { interactions: undefined, height: 200, width: undefined, autoFit: true }) => {
  const c = new Line(dom, {
    data,
    // 默认的x轴依据key为step
    xField,
    // 默认的y轴依据key为data
    yField: 'data',
    // 坐标轴相关
    xAxis: {
      tickCount: 7 // 设置坐标轴刻度数量，防止数据过多导致刻度过密
    },
    yAxis: {
      tickCount: 7,
      min: null
    },
    tooltip: {
      // 如果data为float，则保留4位小数
      formatter: (data) => {
        // console.log(data)
        // FIXME 当前只支持单数据，需要兼容多数据，可以用下面的customContent，但是目前不管
        if (data.data % 1 !== 0) return { name: source[0], value: data.data.toFixed(4) }
        return { name: source[0], value: data.data }
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
    color: defaultColor,
    ...config
  })
  c.render()
  return c
}

// ---------------------------------- 数据格式化 ----------------------------------
/**
 * 为了将数据格式化为图表可用的格式，需要将数据源中的数据进行格式化
 * @param { Object } data 待格式化的数据
 */
const format = (data) => {
  // FIXME 暂时只支持单数据
  const d = data[source[0]].list
  return { d }
}

// ---------------------------------- 渲染、重渲染功能 ----------------------------------
let chartObj = null
// 渲染
const render = (data) => {
  // console.log('渲染折线图')
  const { d } = format(data)
  // console.log('data', data)
  chartObj = createChart(g2Ref.value, d)
}
// 重渲染
const change = (data) => {
  const { d } = format(data)
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
