<template>
  <section
    class="chart-container"
    :class="chartComponent.class"
    @mouseenter="handleMouseEnter"
    @mouseleave="handleMouseLeave"
  >
    <template v-if="status === 'success'">
      <!-- 图表相关控制按钮 -->
      <div class="chart-pannel" v-if="hover && !unknown && !props.chart.error">
        <PannelButton icon="zoom" :tip="$t('experiment.chart.zoom')" @click="zoom" />
      </div>
      <component ref="chartRef" :is="chartComponent.type" :title="chart.name" :chart="chart" />
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
import { onUnmounted, computed, watch } from 'vue'
import { ref } from 'vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
// import SLLoading from '@swanlab-vue/components/SLLoading.vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import LineChart from '@swanlab-vue/charts/LineChart.vue'
import UnknownChart from '@swanlab-vue/charts/UnknownChart.vue'
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
const status = ref(props.chart.error ? 'success' : 'loading')
const unknown = ref(false)
/**
 * @description: 根据chart的type，返回对应的chart组件
 * @return {Object} {type: ChartComponent, class: String}, 其中type为图表组件，class为图表组件的class，class可能不存在
 * @example: {type: LineChart, class: 'line-chart'} 或者 {type: UnknownChart}
 */
const chartComponent = computed(() => {
  switch (props.chart.type) {
    case 'default':
      return {
        type: LineChart
      }
    case 'line':
      return {
        type: LineChart
      }
    default:
      return {
        type: UnknownChart
      }
  }
})

// 监听chartComponent的变化，如果变化为UnknownChart，则unknown为true
// 如果放在chartComponent的computed中会带来额外副作用
watch(
  chartComponent,
  (newVal) => {
    unknown.value = newVal.type === UnknownChart
  },
  {
    immediate: true
  }
)

// ---------------------------------- 订阅 ----------------------------------
let data = {}
// 是否已经渲染，用于控制执行render方法还是change方法
let init = false
const $off = inject('$off')
// 如果props.chart.error存在，则不订阅
props.chart.error ||
  inject('$on')(source, cid, (tag, _tagData, error) => {
    // 异步回调（其实是同步），用于控制图表的显示状态
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
            resolve()
          })
        }
      }
    })
  })

// 卸载时取消订阅
props.chart.error ||
  onUnmounted(() => {
    $off(source, cid)
  })

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
