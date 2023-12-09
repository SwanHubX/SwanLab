<template>
  <div class="chart-container">
    <h1 class="plotly-title absolute z-10">{{ title }}</h1>
    <div class="plotly-chart absolute z-0" ref="gd"></div>
  </div>
</template>

<script setup>
/**
 * @description: 基于 Plotly.js 的图表组件，在此组件中封装
 * @file: PlotlyChart.vue
 * @since: 2023-12-09 00:13:00
 **/

import Plotly from 'plotly.js-dist-min'
import http from '@swanlab-vue/api/http'
import { ref, inject, computed, onUnmounted } from 'vue'

const props = defineProps({
  sources: {
    type: Array,
    required: true
  }
})
// 根据注入的信息，拿到项目id（响应式
const experimentId = inject('experimentId')
// 拿到项目的颜色（响应式
const experimentColor = inject('experimentColor')
// 拿到项目的状态（响应式
const experimentStatus = inject('experimentStatus')

// 根据sources，生成title
const title = computed(() => props.sources.join(' & '))

const gd = ref()

// ---------------------------------- 请求的工具函数 ----------------------------------

/**
 * 根据tag获取到数据
 * @param { string } tag 数据源
 */
const getTag = async (tag) => {
  const { data } = await http.get(`/experiment/${experimentId.value}/` + encodeURIComponent(tag))
  return data.list.map((item) => item.data)
}

// ---------------------------------- 在此处根据sources请求数据 ----------------------------------
;(async function () {
  const data = await getTag(props.sources[0])
  Plotly.newPlot(
    gd.value,
    [{ y: data, line: { color: experimentColor.value } }],
    { showlegend: false },
    { displayModeBar: false, responsive: true }
  )
  // 启动轮询函数，
  startPolling()
})()

// ---------------------------------- 轮询 ----------------------------------
let timer = undefined

const startPolling = () => {
  timer = setInterval(async () => {
    if (experimentStatus.value !== 0) {
      clearInterval(timer)
      return
    }
    const data = await getTag(props.sources[0])
    Plotly.react(gd.value, [{ y: data, line: { color: experimentColor.value } }])
  }, 10000)
}

onUnmounted(() => {
  clearInterval(timer)
})
</script>

<style lang="scss" scoped>
.plotly-chart {
  height: calc(100% + 100px);
  width: calc(100% + 100px);
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
}
.plotly-title {
  top: 10%;
  left: 50%;
  transform: translate(-50%, -50%);
}

.chart-container {
  @apply w-full h-72 border rounded relative overflow-hidden;
  flex: 1 0 400px;
}
</style>
