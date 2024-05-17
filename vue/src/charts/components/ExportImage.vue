<template>
  <SLModal class="pt-5" max-w="890" v-model="downloadModal" escExit>
    <p class="text-lg px-5 font-semibold">
      {{ $t('chart.charts.line.download.export') }} {{ suffix[current].toUpperCase() }}
    </p>
    <div class="flex gap-5 px-5 py-4 md:flex-row flex-col">
      <div class="param">
        <span>{{ $t('chart.charts.line.download.options.name') }}</span>
        <input type="text" class="max-w-56 md:max-w-40" v-model="name" />
      </div>
      <div class="param">
        <span>{{ $t('chart.charts.line.download.options.width') }}</span>
        <input type="number" class="max-w-56 md:max-w-24" v-model="width" />
      </div>
      <div class="param">
        <span>{{ $t('chart.charts.line.download.options.height') }}</span>
        <input type="number" class="max-w-56 md:max-w-24" v-model="height" />
      </div>
    </div>
    <div class="relative p-4" :id="chart.name">
      <LineChart :index="index" ref="downloadRef" :title="chart.name" :chart="chart"></LineChart>
    </div>
    <div class="flex justify-end p-5 pt-0 gap-5">
      <SLMenu menu-class="flex-shrink-0" class="right-0 top-10" down>
        <SLButton class="py-1.5 pl-3 pr-10 rounded">{{ suffix[current].toUpperCase() }}</SLButton>
        <template #pop>
          <SLMenuItems>
            <SLMenuItem v-for="(item, index) in suffix" :key="item" @click="current = index"> {{ item }} </SLMenuItem>
          </SLMenuItems>
        </template>
      </SLMenu>
      <SLButton theme="primary" class="py-1.5 px-3 rounded" @click="download">{{
        $t('chart.charts.line.download.downloadAs', { type: suffix[current].toUpperCase() })
      }}</SLButton>
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
import SLMenu from '@swanlab-vue/components/menu/SLMenu.vue'
import SLMenuItems from '@swanlab-vue/components/menu/SLMenuItems.vue'
import SLMenuItem from '@swanlab-vue/components/menu/SLMenuItem.vue'
import SLButton from '@swanlab-vue/components/SLButton.vue'
import LineChart from '../package/LineChart.vue'
import html2canvas from 'html2canvas'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'

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
const downloadRef = ref(null)

const downloadModal = computed({
  get() {
    downloadRef.value?.render(data)
    props.smoothMethod && downloadRef.value?.smooth(props.smoothMethod)
    addTaskToBrowserMainThread(() => {
      const node = document.getElementById(props.chart.name)
      const legend = node?.querySelector('.lc-legend')
      if (legend) {
        legend.style.overflowY = 'visible'
        legend.style.maxHeight = 'none'
      }
    })
    return props.modelValue
  },
  set(val) {
    emit('update:modelValue', val)
  }
})

const suffix = ['png']
const current = ref(0)
const width = ref(900)
const height = ref(350)
const name = ref(`SwanLab-Chart-${props.chart.name}-${props.index}`)

const download = () => {
  const node = document.getElementById(props.chart.name)
  const legend = node.querySelector('.lc-legend')
  legend?.classList.add('overflow-y-visible')
  html2canvas(node, {
    height: height.value || node.offsetHeight + 10,
    width: width.value || node.offsetWidth + 10,
    scrollY: 0,
    scrollX: 0
    // backgroundColor: 'transparent'
  }).then(async (canvas) => {
    let oImg = new Image()
    oImg.src = canvas.toDataURL()
    const a = document.createElement('a')
    a.download = `${name.value}.${suffix[current.value]}`
    a.href = oImg.src
    a.style.width = `${width.value}`
    a.style.height = `${height.value}`
    a.click()
  })
}
</script>

<style lang="scss" scoped>
.param {
  @apply flex items-center;
  input {
    @apply border px-2 py-1 ml-2 rounded;
  }
  span {
    @apply w-20 md:w-auto block;
  }
}
</style>
