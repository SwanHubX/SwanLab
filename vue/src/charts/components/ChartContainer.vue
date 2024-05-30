<template>
  <section
    class="chart-container"
    :class="chartComponent.class"
    @mouseenter="handleMouseEnter"
    @mouseleave="handleMouseLeave"
  >
    <!-- 错误 -->
    <div class="flex flex-col justify-center grow text-dimmer gap-2" v-if="status === 'error'">
      <SLIcon class="mx-auto h-5 w-5" icon="error" />
      <p class="text-center text-xs">{{ $t('chart.error') }}</p>
    </div>
    <template v-else>
      <!-- 图表相关控制按钮 -->
      <div class="chart-panel">
        <!-- 取消置顶的按钮 -->
        <PanelButton
          class="!text-positive-default"
          icon="pinned"
          click-icon="loading"
          click-class="animate-spin"
          :tip="$t('chart.unpin')"
          @click="unpin"
          v-if="isPinned"
        />
        <template v-if="hover && !unknown && !isAllError">
          <!-- 置顶 -->
          <PanelButton
            icon="pin"
            :tip="$t('chart.pin')"
            @click="pin"
            click-icon="loading"
            click-class="animate-spin"
            v-if="!isPinned"
          />
          <!-- 放大功能 -->
          <PanelButton icon="zoom" :tip="$t('chart.zoom')" @click="zoom" />
          <!-- 更多功能 -->
          <ExportImage :index="index" :chart="chart" :smoothMethod="smoothMethod" v-model="downloadModal" />
          <SLMenu menu-class="flex-shrink-0" class="right-0 top-4">
            <PanelButton icon="more" :tip="$t('chart.more')" />
            <template #pop>
              <SLMenuItems>
                <SLMenuItem class="py-1" @click="unhide" v-if="isHidden">
                  {{ $t('chart.unhide') }}
                </SLMenuItem>
                <SLMenuItem class="py-1" @click="hide" v-else>
                  {{ $t('chart.hide') }}
                </SLMenuItem>
                <SLMenuItem class="py-1" v-if="isLineChart" @click="downloadModal = true">
                  {{ $t('chart.charts.share.download.button') }}
                </SLMenuItem>
              </SLMenuItems>
            </template>
          </SLMenu>
        </template>
      </div>
      <component :index="index" ref="chartRef" :is="chartComponent.type" :title="chart.name" :chart="chart" />
    </template>
    <!-- 加载中 -->
    <div class="h-full w-full top-0 left-0 bg-default flex items-center justify-center absolute rounded" v-if="loading">
      <SLLoading />
    </div>
  </section>
</template>

<script setup>
/**
 * @description: 图表容器，用于包裹并渲染图表组件
 * @file: ChartContainer.vue
 * @since: 2023-12-12 21:01:26
 **/
import { inject } from 'vue'
import { onUnmounted, computed, watch, onMounted } from 'vue'
import { ref } from 'vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
// import SLLoading from '@swanlab-vue/components/SLLoading.vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import LineChart from '../package/LineChart.vue'
import AudioChart from '../package/AudioChart.vue'
import TextChart from '../package/TextChart.vue'
import UnknownChart from '../package/UnknownChart.vue'
import PanelButton from './PanelButton.vue'
import { debounce } from '@swanlab-vue/utils/common'
import ImageChart from '../package/ImageChart.vue'
import SLLoading from '@swanlab-vue/components/SLLoading.vue'
import SLMenu from '@swanlab-vue/components/menu/SLMenu.vue'
import SLMenuItems from '@swanlab-vue/components/menu/SLMenuItems.vue'
import SLMenuItem from '@swanlab-vue/components/menu/SLMenuItem.vue'
import ExportImage from '../components/ExportImage.vue'

const props = defineProps({
  chart: {
    type: Object,
    required: true
  },
  index: {
    type: Number,
    required: true
  },
  isPinned: {
    type: Boolean,
    required: true
  },
  isHidden: {
    type: Boolean,
    required: true
  }
})

const cid = props.chart.id
const source = props.chart.source
const emit = defineEmits(['pin', 'unpin', 'hide', 'unhide'])

// ---------------------------------- 判断是否处于hover状态 ----------------------------------
const hover = ref(false)
const handleMouseEnter = () => (hover.value = true)
const handleMouseLeave = () => (hover.value = false)

// ---------------------------------- 控制chart显示状态 ----------------------------------
// 是否全部数据都错
const isAllError = computed(() => {
  return Object.keys(props.chart.error || {}).length === source.length
})

const status = ref(props.chart.error ? 'success' : undefined)
const loading = ref(!isAllError.value)

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
        type: LineChart,
        class: 'line-chart'
      }
    case 'line':
      return {
        type: LineChart,
        class: 'line-chart'
      }
    case 'audio':
      return {
        type: AudioChart,
        class: 'audio-chart'
      }
    case 'image':
      return {
        type: ImageChart,
        class: 'image-chart'
      }
    case 'text':
      return {
        type: TextChart,
        class: 'text-chart'
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

// 判断是否为Promise或者AsyncFunction
function isPromiseAndAsyncFunction(func) {
  return (
    (func instanceof Promise ||
      (func !== null &&
        func !== undefined &&
        typeof func === 'object' &&
        typeof func.then === 'function' &&
        typeof func.catch === 'function')) &&
    func.constructor &&
    func.constructor.name === 'AsyncFunction'
  )
}
// ---------------------------------- 订阅 ----------------------------------
let data = {}
provide('data', data)
// 是否已经渲染，用于控制执行render方法还是change方法
let hasInited = false
const $off = inject('$off')
// 如果props.chart.error存在，则不订阅
onMounted(() => {
  // error的key数量不等于source的长度
  !isAllError.value &&
    inject('$on')(source, cid, (tag, _tagData, error) => {
      // 异步回调（其实是同步），用于控制图表的显示状态
      return new Promise((resolve, reject) => {
        if (error) {
          status.value = 'error'
          loading.value = false
          reject(error)
        } else {
          data[tag] = _tagData
          // 判断data的key数量是否和source长度相同
          if (Object.keys(data).length === source.length) {
            // 渲染
            addTaskToBrowserMainThread(() => {
              if (!hasInited) {
                // 判断render方法是否为Promise，如果是，则等待渲染完成再resolve并设置status为success
                render(data)
              } else change(data)
              hasInited = true
              resolve()
            })
          }
        }
      })
    })
})

const isChartError = computed(() => {
  return props.chart.error && Object.keys(props.chart.error).length > 0
})

// 卸载时取消订阅
isChartError.value ||
  onUnmounted(() => {
    $off(source, cid)
  })

// ---------------------------------- 图表对象控制 ----------------------------------
const chartRef = ref(null)
// 渲染功能
const render = debounce(() => {
  // 判断是否为Promise，如果是，则等待渲染完成再resolve并设置loading为false
  // console.log('isPromise(chartRef.value.render)', isPromise(chartRef.value.render), chartRef.value.render)
  const r = chartRef.value?.render
  if (isPromiseAndAsyncFunction(r)) r(data).finally(() => (loading.value = false))
  else {
    r && r(data)
    loading.value = false
  }
}, 100)
// 重渲染功能
const change = debounce(() => {
  chartRef.value.change(data)
}, 100)
// 放大功能
const zoom = () => {
  chartRef.value.zoom(data)
}

// ---------------------------------- 图表显示/隐藏tooltip ----------------------------------

const showTooltip = (...args) => {
  chartRef.value?.showTooltip && chartRef.value.showTooltip(...args)
}

const hideTooltip = () => {
  chartRef.value?.hideTooltip && chartRef.value.hideTooltip()
}

// ---------------------------------- 图表加粗 ----------------------------------
const thicken = (...args) => {
  chartRef.value?.thicken && chartRef.value.thicken(...args)
}

const thin = (...args) => {
  chartRef.value?.thin && chartRef.value.thin(...args)
}

// ---------------------------------- 图表平滑方法 ----------------------------------
const smoothMethod = ref(null)
const smooth = (method) => {
  smoothMethod.value = method
  chartRef.value?.smooth && chartRef.value.smooth(method)
}

// ---------------------------------- 图表置顶 ----------------------------------
const pin = () => {
  emit('pin', props.chart)
}
const unpin = () => {
  emit('unpin', props.chart)
}

// ---------------------------------- 图表隐藏 ----------------------------------

const hide = () => {
  emit('hide', props.chart)
}

const unhide = () => {
  emit('unhide', props.chart)
}

// ---------------------------------- 下载折线图 ----------------------------------

const downloadModal = ref(false)

const isLineChart = computed(() => {
  if (isChartError.value) return false
  return props.chart?.type === 'default' || props.chart?.type === 'line'
})

// ---------------------------------- 暴露组件对象 ----------------------------------
defineExpose({
  chartRef,
  smooth,
  thicken,
  thin,
  showTooltip,
  hideTooltip
})
</script>

<style lang="scss" scoped>
.chart-container {
  @apply w-full h-auto min-h-[288px] border rounded relative bg-default;
  @apply px-3 py-4;
  @apply flex-col flex justify-between;
}

.chart-panel {
  @apply absolute top-1 right-2 h-4;
  @apply flex justify-end items-center gap-2;
}
// ---------------------------------- 图表样式 ----------------------------------
.line-chart {
  @apply col-span-3 lg:col-span-1;
}
.audio-chart,
.image-chart,
.text-chart {
  @apply col-span-3;
}
</style>
