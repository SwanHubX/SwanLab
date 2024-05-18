<template>
  <SLModal max-w="890" class="overflow-hidden" v-model="downloadModal" escExit>
    <h1 class="text-lg px-5 font-semibold py-4 border-b">
      {{ $t('chart.charts.share.download.png.title') }}
    </h1>
    <div class="flex gap-5 px-5 py-6 md:flex-row flex-col">
      <div class="param">
        <span>{{ $t('chart.charts.share.download.png.name') }}:</span>
        <input
          type="text"
          class="w-40 text-sm border px-2 py-1 rounded ml-2"
          :placeholder="placeholder"
          v-model="name"
        />
      </div>
      <div class="param">
        <span>{{ $t('chart.charts.share.download.png.width') }}:</span>
        <InputSlideBar min="200" max="1200" v-model="width" />
      </div>
      <div class="param">
        <span>{{ $t('chart.charts.share.download.png.height') }}:</span>
        <InputSlideBar min="200" max="1200" v-model="height" />
      </div>
    </div>
    <!-- 图像主体 -->
    <div class="mx-auto mb-6 py-4" :style="{ width: `${width}px` }">
      <div class="relative w-full h-full" ref="chartRef">
        <LineChart :index="index" ref="lineChartRef" :title="chart.name" :chart="chart" />
      </div>
    </div>
    <div class="flex justify-end p-5 pt-0 gap-5">
      <SLButton
        theme="primary"
        class="py-1.5 px-3 rounded"
        @click="download"
        :text="$t('chart.charts.share.download.png.download')"
      />
    </div>
  </SLModal>
</template>

<script setup>
/**
 * @description:
 * @file: ExportImage.vue
 * @since: 2024-05-14 16:22:09
 **/
import SLModal from '@swanlab-vue/components/SLModal.vue'
import SLButton from '@swanlab-vue/components/SLButton.vue'
import LineChart from '../package/LineChart.vue'
import { generateFileNameWithTime } from '@swanlab-vue/utils/common'
import html2canvas from 'html2canvas'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import InputSlideBar from './InputSlideBar.vue'

const props = defineProps({
  chart: {
    type: Object,
    required: true
  },
  index: {
    type: Number,
    required: true
  },
  modelValue: {
    type: Boolean,
    default: false
  },
  smoothMethod: {
    type: Object,
    default: null
  }
})
const emit = defineEmits(['update:modelValue'])
const data = inject('data')
/**
 * @type {import('vue').Ref} lineChartRef
 */
const lineChartRef = ref(null)
/**
 * @type {import('vue').Ref<HTMLDivElement>} chartRef
 */
const chartRef = ref(null)

const downloadModal = computed({
  get() {
    // 显示时渲染图表
    if (props.modelValue) {
      lineChartRef.value?.render(data)
      props.smoothMethod && lineChartRef.value?.smooth(props.smoothMethod)
    } else {
      // 关闭时清空长宽数值
      width.value = originalWidth
      height.value = originalHeight
    }
    return props.modelValue
  },
  set(val) {
    emit('update:modelValue', val)
  }
})

// ---------------------------------- 名称、长宽等属性 ----------------------------------

const placeholder = generateFileNameWithTime('SwanLab-Chart').split('.')[0]
const name = ref('')
const originalWidth = 600
const originalHeight = 400
const width = ref(originalWidth)
const height = ref(originalHeight)

// ---------------------------------- 设置动态拦截高度设置 ----------------------------------

watch([height, downloadModal, lineChartRef], (newVal) => {
  // 设置高度
  if (!lineChartRef.value) return
  if (!newVal[1]) return
  const h = newVal[0]
  console.log(h)
  /**
   * @type {setOriginalChartHeight} lineChartRef.value.setHeight
   */
  addTaskToBrowserMainThread(() => {
    lineChartRef.value.setHeight(h)
  })
})

// ---------------------------------- 下载按钮 ----------------------------------

const download = () => {
  html2canvas(chartRef.value).then(async (canvas) => {
    // 新建canvas，为canvas本身增加一圈白色边框
    const borderWidth = 32
    const newCanvas = document.createElement('canvas')
    const ctx = newCanvas.getContext('2d')
    newCanvas.width = canvas.width + borderWidth * 2
    newCanvas.height = canvas.height + borderWidth * 2
    ctx.fillStyle = '#fff'
    ctx.fillRect(0, 0, newCanvas.width, newCanvas.height)
    ctx.drawImage(canvas, borderWidth, borderWidth)
    // 创建一个临时链接元素
    const link = document.createElement('a')
    link.href = newCanvas.toDataURL('image/png')
    link.download = `${name.value || placeholder}.png`
    link.click()
  })
}
</script>

<style lang="scss" scoped>
.param {
  @apply flex items-center text-dimmer;
  span {
    @apply w-20 md:w-auto block;
  }
}
</style>
