<template>
  <div class="chart-container">
    <p class="text-center">{{ title }}</p>
    <div ref="g2"></div>
  </div>
</template>

<script setup>
/**
 * @description: 基于@antv/g2plot的图表组件
 * @file: G2Chart.vue
 * @since: 2023-12-11 10:45:40
 **/
import { Line } from '@antv/g2plot'
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

const g2 = ref()

// ---------------------------------- 请求的工具函数 ----------------------------------

/**
 * 根据tag获取到数据
 * @param { string } tag 数据源
 */
const getTag = async (tag) => {
  const { data } = await http.get(`/experiment/${experimentId.value}/` + encodeURIComponent(tag))
  // FIXME 为data添加一个step字段
  data.list.forEach((item, index) => {
    item.step = index
  })
  return data.list
}

// ---------------------------------- 在此处根据sources请求数据 ----------------------------------
let chart
;(async function () {
  const data = await getTag(props.sources[0])
  chart = new Line(g2.value, {
    data,
    xField: 'step',
    yField: 'data',
    // 坐标轴相关
    xAxis: {
      // type: 'timeCat',
      tickCount: 5 // 设置坐标轴刻度数量，防止数据过多导致刻度过密
    },
    yAxis: {},
    // 大小相关
    height: 200,
    autoFit: true,
    // 开启一些交互
    interactions: [{ type: 'brush-x' }],
    // 样式相关
    // smooth: true, // 平滑曲线
    color: experimentColor.value
  })

  chart.render()
  // 启动轮询函数，
  startPolling()
})()

// ---------------------------------- 轮询 ----------------------------------
let timer = undefined

const startPolling = () => {
  timer = setInterval(async () => {
    if (experimentStatus.value !== 0) return clearInterval(timer)
    const data = await getTag(props.sources[0])
    chart.changeData(data)
  }, 1000)
}

onUnmounted(() => {
  clearInterval(timer)
})
</script>

<style lang="scss" scoped>
.chart-container {
  @apply w-full h-72 border rounded relative overflow-hidden;
  @apply px-3 py-4;
  @apply flex-col flex justify-between;
  flex: 1 0 400px;
}
</style>
