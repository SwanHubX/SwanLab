<template>
  <ChartContainer :title="title">
    <template #pannel>
      <PannelButton icon="zoom" :tip="$t('experiment.chart.zoom')" @click="handleOpen" />
    </template>
    <!-- x轴坐标单位 -->
    <p class="absolute right-5 bottom-10 text-xs text-dimmer scale-90">step</p>
    <div ref="g2Ref"></div>
    <!-- 放大效果 -->
    <SLModal class="p-10 pt-0" max-w="-1" v-model="zoom">
      <p class="text-center mt-4 mb-10 text-2xl font-semibold">{{ title }}</p>
      <div ref="g2ZoomRef"></div>
    </SLModal>
  </ChartContainer>
</template>

<script setup>
/**
 * @description: 基于@antv/g2plot的图表组件
 * @file: G2Chart.vue
 * @since: 2023-12-11 10:45:40
 **/
import { Line } from '@antv/g2plot'
import http from '@swanlab-vue/api/http'
import ChartContainer from './ChartContainer.vue'
import { ref, inject, computed, onUnmounted } from 'vue'
import PannelButton from './PannelButton.vue'
import SLModal from '@swanlab-vue/components/SLModal.vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import { useExperimentStroe } from '@swanlab-vue/store'
const experimentStore = useExperimentStroe()
const props = defineProps({
  sources: {
    type: Array,
    required: true
  }
})
// 项目id
const id = computed(() => experimentStore.id)
// 拿到图表颜色
const color = experimentStore.defaultColor
// 拿到项目的状态（响应式
const status = computed(() => experimentStore.status)

// 根据sources，生成title
const title = computed(() => props.sources.join(' & '))

const g2Ref = ref()

// ---------------------------------- 请求的工具函数 ----------------------------------
const tagData = ref()
/**
 * 根据tag获取到数据
 * @param { string } tag 数据源
 */
const getTag = async (tag) => {
  const { data } = await http.get(`/experiment/${id.value}/tag/` + encodeURIComponent(tag))
  // FIXME 为data添加一个step字段
  data.list.forEach((item, index) => {
    item.step = index
  })
  tagData.value = data
  return data.list
}

// ---------------------------------- 在此处根据sources请求数据 ----------------------------------
let chart
const createChart = (dom, data, interactions = undefined, height = 200, width = undefined, autoFit = true) => {
  const c = new Line(dom, {
    data,
    xField: 'step',
    yField: 'data',
    // 坐标轴相关
    xAxis: {
      // type: 'timeCat',
      text: 'x 轴标题',
      tickCount: 7 // 设置坐标轴刻度数量，防止数据过多导致刻度过密
    },
    yAxis: {
      tickCount: 7
    },
    // 大小相关
    height,
    width,
    autoFit,
    // 开启一些交互
    interactions,
    // 样式相关
    // smooth: true, // 平滑曲线
    color
  })
  c.render()
  return c
}
;(async function () {
  const data = await getTag(props.sources[0])
  chart = createChart(g2Ref.value, data, undefined)
  // 启动轮询函数，
  startPolling()
})()

// ---------------------------------- 轮询 ----------------------------------
let timer = undefined

const startPolling = () => {
  timer = setInterval(async () => {
    if (status.value !== 0) return clearInterval(timer)
    const data = await getTag(props.sources[0])
    chart.changeData(data)
  }, 1000)
}

onUnmounted(() => {
  clearInterval(timer)
})

// ---------------------------------- 控制放大 ----------------------------------
const zoom = ref(false)
const g2ZoomRef = ref()
const handleOpen = () => {
  zoom.value = true
  // 当前window的高度
  const height = window.innerHeight * 0.6
  addTaskToBrowserMainThread(() => {
    createChart(
      g2ZoomRef.value,
      tagData.value.list,
      [{ type: 'view-zoom' }, { type: 'drag-move' }],
      height,
      undefined,
      true
    )
  })
}
</script>

<style lang="scss" scoped></style>
