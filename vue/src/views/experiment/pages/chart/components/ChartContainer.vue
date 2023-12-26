<template>
  <section class="chart-container" @mouseenter="handleMouseEnter" @mouseleave="handleMouseLeave">
    <template v-if="status === 'success'">
      <!-- 图表相关控制按钮 -->
      <div class="chart-pannel" v-if="hover">
        <PannelButton icon="zoom" :tip="$t('experiment.chart.zoom')" @click="zoom" />
      </div>
      <component ref="chartRef" :is="chartComponent(chart.type)" :title="title(chart.source[0])" :chart="chart" />
    </template>
    <!-- 错误 -->
    <div class="flex flex-col justify-center grow text-dimmer gap-2" v-else-if="status === 'error'">
      <SLIcon class="mx-auto h-5 w-5" icon="error" />
      <p class="text-center text-xs">{{ $t('experiment.chart.error') }}</p>
    </div>
    <!-- 加载中，暂无加载动画 -->
    <!-- <div class="h-full flex items-center justify-center" v-else>
      <SLLoading />
    </div> -->
  </section>
</template>

<script setup>
/**
 * @description: 图表容器，用于包裹并渲染图表组件
 * @file: ChartContainer.vue
 * @since: 2023-12-12 21:01:26
 **/
import { inject } from 'vue'
import { onUnmounted } from 'vue'
import { ref } from 'vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
// import SLLoading from '@swanlab-vue/components/SLLoading.vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import LineChart from './package/LineChart.vue'
import PannelButton from './PannelButton.vue'
import { debounce } from '@swanlab-vue/utils/common'
const props = defineProps({
  chart: {
    type: Object,
    required: true
  }
})

const cid = props.chart._cid
const source = props.chart.source

// ---------------------------------- 判断是否处于hover状态 ----------------------------------
const hover = ref(false)
const handleMouseEnter = () => (hover.value = true)
const handleMouseLeave = () => (hover.value = false)

// ---------------------------------- 控制chart显示状态 ----------------------------------
const status = ref('loading')

const chartComponent = (type) => {
  switch (type) {
    case 'default':
      return LineChart
  }
}

// ---------------------------------- 订阅 ----------------------------------
let data = {}
// 是否已经渲染
let init = false
const $off = inject('$off')
inject('$on')(source, cid, (tag, _tagData, error) => {
  return new Promise((resolve, reject) => {
    if (error) {
      status.value = 'error'
      reject(error)
    } else {
      status.value = 'success'
      data[tag] = _tagData
      // 判断data的key数量是否和source长度相同
      if (Object.keys(data).length === source.length) {
        // 渲染
        addTaskToBrowserMainThread(() => {
          if (!init) render(data)
          else change(data)
          init = true
        })
      }
      resolve(data)
    }
  })
})

// 卸载时取消订阅
onUnmounted(() => {
  $off(source, cid)
})

// ---------------------------------- 图表标题控制 ----------------------------------
const title = (t) => {
  // TODO 多数据时可能有不同的显示方式
  return t
}

// ---------------------------------- 图表对象控制 ----------------------------------
const chartRef = ref(null)
// 渲染功能
const render = debounce(() => {
  chartRef.value.render(data)
}, 100)
// 重渲染功能
const change = debounce(() => {
  chartRef.value.change(data)
}, 100)

// 放大功能
const zoom = () => {
  chartRef.value.zoom(data)
}
</script>

<style lang="scss" scoped>
.chart-container {
  @apply w-full h-72 border rounded relative overflow-hidden;
  @apply px-3 py-4;
  @apply flex-col flex justify-between;
}

.chart-pannel {
  @apply absolute top-1 right-2   h-4;
  @apply flex justify-end gap-2;
}
</style>
