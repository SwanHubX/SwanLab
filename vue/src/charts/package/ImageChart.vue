<template>
  <!-- 图表标题 -->
  <p class="text-center font-semibold">{{ title }}</p>
  <!-- 如果图表数据错误 -->
  <div class="flex flex-col justify-center grow text-dimmer gap-2" v-if="error">
    <SLIcon class="mx-auto h-5 w-5" icon="error" />
    <p class="text-center text-xs">
      <!-- 在此处显示错误信息 -->
      {{ $t('common.chart.charts.image.error', { type: error['data_class'], tag: source[0] }) }}
    </p>
  </div>
  <!-- 如果图表数据正确 -->
  <template v-else>
    <!-- 在此处完成图表主体定义 -->
    <div class="image-content image-content-no-zoom" ref="imageContentRef">
      <div class="flex flex-col justify-center items-center h-full" v-if="loading">
        <SLLoading />
      </div>
      <!-- 加载完成 -->
      <div class="images-container" :style="setGrid(stepsData[currentIndex][source[0]].length)" v-else>
        <div
          class="flex flex-col items-center justify-center relative"
          v-for="(s, index) in stepsData[currentIndex][source[0]]"
          :key="index"
        >
          <div class="image-container">
            <img :src="imagesData[s.filename].url" @click="handelClickZoom(s.filename)" />
            <DownloadButton class="!w-5 !h-5 !p-0.5 download-button" @click.stop="download(s.filename)" />
          </div>
          <p class="text-xs">{{ s.caption }}</p>
        </div>
      </div>
    </div>
    <div class="h-8">
      <SlideBar
        class="mt-2"
        v-model="currentIndex"
        :max="maxIndex"
        :min="minIndex"
        :bar-color="barColor"
        :key="slideKey"
        v-if="maxIndex !== minIndex"
      />
    </div>
    <!-- 放大效果弹窗 -->
    <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom">
      <p class="text-center mt-4 mb-10 text-2xl font-semibold">{{ title }}</p>
      <div class="image-content flex flex-col justify-center">
        <div class="flex flex-col justify-center items-center" v-if="loading">
          <SLLoading />
        </div>
        <!-- 加载完成 -->
        <div class="images-container" :style="setGrid(stepsData[currentIndex][source[0]].length)" v-else>
          <div
            class="flex flex-col items-center h-full justify-center"
            v-for="(s, index) in stepsData[currentIndex][source[0]]"
            :key="index"
          >
            <div class="image-container">
              <img :src="imagesData[s.filename].url" @click="handelClickZoom(s.filename)" />
              <DownloadButton class="!w-5 !h-5 !p-0.5 download-button" @click.stop="download(s.filename)" />
            </div>
            <p class="text-xs mt-2">{{ s.caption }}</p>
          </div>
        </div>
      </div>
      <div class="h-8">
        <SlideBar
          class="mt-2"
          v-model="currentIndex"
          :max="maxIndex"
          :min="minIndex"
          :bar-color="barColor"
          :key="slideKey"
          v-if="maxIndex !== minIndex"
        />
      </div>
    </SLModal>
    <!-- 额外的放大功能，点击某个图像，放大显示 -->
    <SLModal
      class="w-full flex justify-center min-h-[calc(100vh-8rem)] p-14"
      max-w="-1"
      v-model="isSingleZoom"
      close-on-overlay-click
    >
      <img :src="imagesData[signleZoomFilename].url" class="object-contain" />
    </SLModal>
  </template>
</template>

<script setup>
/**
 * @description: 图像图表组件
 * @file: ImageChart.vue
 * @since: 2024-02-03 13:24:57
 **/
import SLModal from '@swanlab-vue/components/SLModal.vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import { ref, inject, watch, computed } from 'vue'
import { useExperimentStore } from '@swanlab-vue/store'
import SlideBar from '../components/SlideBar.vue'
import * as UTILS from './utils'
import { debounce } from '@swanlab-vue/utils/common'
import DownloadButton from '../components/DownloadButton.vue'
const experimentStore = useExperimentStore()
// ---------------------------------- 配置 ----------------------------------
const props = defineProps({
  title: {
    type: String,
    required: true
  },
  chart: {
    type: Object,
    required: true
  },
  index: {
    type: Number,
    required: true
  }
})
const run_id = computed(() => experimentStore.experiment.run_id)
const source = computed(() => {
  return props.chart.source
})

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------
// TODO 当前只支持单个tag，所以error就是error.{tag}
const error = ref(props.chart.error[source.value[0]])

// ---------------------------------- 图表颜色配置 ----------------------------------
// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const colors = inject('colors')
if (!colors) throw new Error('colors is not defined, please provide colors in parent component')
// ---------------------------------- 组件渲染逻辑 ----------------------------------
// 已经滑动部分颜色，应该通过色盘计算得到
const barColor = inject('colors')[0]
// 当前滑块索引
const __currentIndex = ref(0)
// 最小索引
const minIndex = ref(undefined)
// 最大索引
const maxIndex = ref(undefined)
// slide的key
const slideKey = ref(0)
watch([maxIndex, minIndex], ([max, min]) => {
  slideKey.value = max + '-' + min
})
const loading = ref(true)
const imageContentRef = ref(null)

// 计算属性，拦截滑块索引变化
const currentIndex = computed({
  get: () => __currentIndex.value,
  set: (val) => {
    // 如果当前值存在于stepsData的key中，则直接赋值
    if (val in stepsData) {
      if (val === __currentIndex.value && loading.value === false) return
      __currentIndex.value = val
    } else {
      // 寻找一个最接近的值
      const keys = Array.from(Object.keys(stepsData))
      const index = keys.findIndex((item) => item > val)
      const num = Number(index === -1 ? keys[keys.length - 1] : keys[index])
      if (num === __currentIndex.value) return
      __currentIndex.value = num
    }
    loading.value = true
    imageContentRef.value.height = imageContentRef.value.offsetHeight + 'px'
    debounceGetImagesData(stepsData[__currentIndex.value])
  }
})

// 布局处理,一共length列，最多显示8列
const setGrid = (length) => {
  if (length <= 8) {
    return 'grid-template-columns:' + 'repeat(' + length + ', minmax(0, 1fr))'
  } else {
    return 'grid-template-columns:' + 'repeat(8, minmax(0, 1fr))'
  }
}

// ---------------------------------- 数据格式化 ----------------------------------
// 前端映射数据，不包含图像，数据格式：{step: {tag: [{filename: string, caption: string}]}}
const stepsData = {}
// 图像数据缓存, 数据格式：{String<filename>: {blob: Blob, url: string<base64>}}
const imagesData = {}

const getImagesData = async (stepData) => {
  // 单数据
  const tag = source.value[0]
  const promises = []
  // console.log('stepData', stepData)
  for (const { filename } of stepData[tag]) {
    if (!imagesData[filename]) {
      promises.push(
        new Promise((resolve) => {
          UTILS.media.get(filename, run_id.value, tag).then((blob) => {
            // blob转换为图像base64
            const reader = new FileReader()
            reader.onload = function () {
              imagesData[filename] = { blob, url: reader.result }
              resolve()
            }
            reader.readAsDataURL(blob)
          })
        })
      )
    }
  }
  await Promise.all(promises)
  loading.value = false
  imageContentRef.value.height = ''
  // console.log('图像数据：', imagesData)
}

const debounceGetImagesData = debounce(getImagesData, 500)

/**
 * 将后端传回的数据转换为图像数据，存储到imageStepsData中
 * 后端返回的数据格式为：{tag: {..., list:[{data: string or list<strinf>, more: object or list<object>, index: string<number>, create_time: string}]}}
 * @param { Object } data 后端返回的数据
 */
const changeData2Image = (data) => {
  // 遍历数据，将数据转换为图像数据
  let _maxIndex = 0
  let _minIndex = Infinity
  for (const tag in data) {
    // 遍历tag下的数据
    _maxIndex = Math.max(Number(data[tag].list[data[tag].list.length - 1].index), _maxIndex)
    _minIndex = Math.min(Number(data[tag].list[0].index), _minIndex)
    for (const item of data[tag].list) {
      if (!stepsData[item.index]) stepsData[item.index] = {}
      else continue
      // 对当前tag下的数据进行处理,向stepsData中添加数据
      if (!stepsData[item.index][tag]) stepsData[item.index][tag] = []
      // 添加数据,如果data是字符串，则直接添加，如果是数组，则遍历添加
      if (typeof item.data === 'string') {
        stepsData[item.index][tag].push({ filename: item.data, caption: item.more?.caption })
      } else {
        for (let i = 0; i < item.data.length; i++) {
          stepsData[item.index][tag].push({ filename: item.data[i], caption: item.more[i]?.caption })
        }
      }
    }
  }
  maxIndex.value = _maxIndex
  minIndex.value = _minIndex
}

// ---------------------------------- 渲染、重渲染功能 ----------------------------------

// 渲染
const render = (data) => {
  // 数据格式化
  changeData2Image(data)
  // console.log('图像数据：', stepsData)
  currentIndex.value = minIndex.value
}
// 重渲染
const change = (data) => {
  // 数据格式化
  changeData2Image(data)
}

// ---------------------------------- 放大功能 ----------------------------------
// 是否放大
const isZoom = ref(false)
// 放大数据
const zoom = () => {
  isZoom.value = true
}
// 放大时，左右键切换
const handelKeydown = (e) => {
  // 如果点了esc，则关闭放大
  if (e.key === 'Escape') {
    isZoom.value = false
    return
  }
  if (e.key !== 'ArrowLeft' && e.key !== 'ArrowRight') return
  // console.log('e', e)
  // console.log('key', e.key)
  const steps = Object.keys(stepsData)
  const index = steps.findIndex((item) => item == currentIndex.value)
  if (e.key === 'ArrowLeft') {
    // 向左
    // console.log('减小', stepsData)
    if (currentIndex.value == minIndex.value) return
    currentIndex.value = Number(steps[index - 1])
  } else {
    // 向右
    // console.log('增大', stepsData)
    if (currentIndex.value == maxIndex.value) return
    currentIndex.value = Number(steps[index + 1])
  }
  // console.log('currentIndex', currentIndex.value)
}

watch(isZoom, (val) => {
  if (val) {
    document.addEventListener('keydown', handelKeydown)
  } else {
    document.removeEventListener('keydown', handelKeydown)
  }
})
// ---------------------------------- 点击某个图像，放大 ----------------------------------
const isSingleZoom = ref(false)
watch(isSingleZoom, () => {
  isZoom.value = false
})
const signleZoomFilename = ref()
// 点击某个图像，放大
const handelClickZoom = (filename) => {
  signleZoomFilename.value = filename
  isSingleZoom.value = true
  // console.log('filename', filename)
  // console.log('image data', imagesData[filename])
}

// ---------------------------------- 点击下载按钮下载 ----------------------------------

const download = (filename) => {
  // 拿到blob数据
  const blob = imagesData[filename].blob
  // 创建a标签
  const a = document.createElement('a')
  // 创建url
  const url = window.URL.createObjectURL(blob)
  // 设置a标签
  a.href = url
  a.download = filename
  // 模拟点击
  a.click()
}

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  render,
  change,
  zoom
})
</script>

<style lang="scss" scoped>
.image-content {
  @apply mt-1 p-2 w-full rounded-sm relative min-h-[224px];
  img {
    @apply cursor-pointer;
  }
  .images-container {
    @apply grid gap-2 h-full;
    .image-container {
      @apply inline-block relative;

      &:hover .download-button {
        @apply block;
      }
      .download-button {
        @apply opacity-75 absolute top-1 right-1 hidden;
      }
    }
  }
}

.image-content-no-zoom {
  @apply h-56 overflow-y-auto overflow-x-clip;
}
</style>
