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
      <!-- 单实验图表 -->
      <div
        class="images-container"
        :style="setGrid(stepsData[currentIndex][source[0]].length)"
        v-if="!loading && !isMulti"
      >
        <div class="image-detail" v-for="(s, index) in stepsData[currentIndex][source[0]]" :key="index">
          <div class="image-container">
            <img :src="imagesData[s.filename].url" @click="handelClickZoom(s.filename)" />
            <DownloadButton class="download-button" @click.stop="download(s.filename)" />
          </div>
          <p class="text-xs">{{ s.caption }}</p>
        </div>
      </div>
      <!-- 多实验图表 -->
      <div v-if="!loading && isMulti" class="images-container" :style="setGrid(visiableSources.length)">
        <div class="image-detail" v-for="(s, name) in stepsData[currentIndex]" :key="name">
          <div class="text-xs text-left w-full flex items-center pb-1">
            <div class="h-2 w-2 rounded-full border"></div>
            <p class="pl-2">{{ name }}</p>
          </div>
          <div class="image-container">
            <img
              :src="imagesData[s[currentInnerIndex - 1].filename].url"
              @click="handelClickZoom(s[currentInnerIndex - 1].filename)"
            />
            <DownloadButton class="download-button" @click.stop="download(s[currentInnerIndex - 1].filename)" />
          </div>
          <p class="text-xs">{{ s.caption }}</p>
        </div>
      </div>
    </div>
    <div class="h-8 flex items-center justify-center">
      <SlideBar
        class="mt-2"
        v-model="currentIndex"
        :max="maxIndex"
        :min="minIndex"
        :bar-color="barColor"
        :key="slideKey"
        @turn="handelTurn"
        v-if="maxIndex !== minIndex"
      />
      <SlideBar
        v-model="currentInnerIndex"
        :max="maxInnerIndex"
        :min="minInnerIndex"
        :bar-color="barColor"
        :key="maxInnerIndex"
        reference="index"
        v-if="source?.length > 1"
      />
    </div>
    <!-- 放大效果弹窗 -->
    <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom">
      <p class="text-center mt-4 mb-10 text-2xl font-semibold">{{ title }}</p>
      <div class="image-content image-content-zoom">
        <div class="flex flex-col justify-center items-center h-full" v-if="loading">
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
              <DownloadButton class="download-button" @click.stop="download(s.filename)" />
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
          @turn="handelTurn"
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
      <div class="image-single-zoom">
        <div class="relative">
          <img :src="imagesData[signleZoomFilename].url" class="object-contain" />
          <DownloadButton class="download-button" @click.stop="download(signleZoomFilename)" />
        </div>
      </div>
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
import { useExperimentStore, useProjectStore } from '@swanlab-vue/store'
import SlideBar from '../components/SlideBar.vue'
import * as UTILS from './utils'
import { debounce } from '@swanlab-vue/utils/common'
import DownloadButton from '../components/DownloadButton.vue'
import { useRoute } from 'vue-router'

const experimentStore = useExperimentStore()
const projectStore = useProjectStore()
const route = useRoute()

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
/**
 * 图表数据来源，为数组
 * 当isMulti为true时，数组元素为exp，即给该表提供数据的实验的名称
 * 当isMulti为false时，数组元素为tag，即给该表提供数据的tag
 */
const source = computed(() => {
  return props.chart.source
})
// 是否为多实验的图表，根据路由名称判断
const isMulti = computed(() => {
  return route.name === 'charts'
})
/**
 * 实验名对应的run_id
 * 单实验时直接从 experiment pinia 中拿
 * 多实验时通过 project pinia 中的 getExpRunIdByName 获取（传入实验名查询 run_id）
 */
const run_id = computed(() => {
  if (!isMulti.value) return experimentStore.experiment?.run_id
  const run_ids = {}
  for (const exp of source.value) {
    run_ids[exp] = projectStore.getExpRunIdByName(exp)
  }
  return run_ids
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
    // 如果是多实验模式，将内部数据索引置 1
    if (isMulti.value) currentInnerIndex.value = 1
    // 获取数据
    debounceGetImagesData(stepsData[__currentIndex.value])
    console.log(stepsData)
  }
})

// 事件处理，触发slideBar的turn事件
const handelTurn = (direction, value) => {
  const keys = Array.from(Object.keys(stepsData))
  const index = keys.findIndex((item) => item > value)
  if (direction === 'forward') {
    currentIndex.value = Number(index === -1 ? keys[keys.length - 1] : keys[index])
  } else {
    // 向下获取
    if (index === -1) currentIndex.value = Number(keys[Math.max(0, keys.length - 2)])
    else currentIndex.value = Number(keys[Math.max(0, index - 2)])
  }
}

// 布局处理,一共length列，最多显示8列
const setGrid = (length) => {
  if (length <= 8) {
    return 'grid-template-columns:' + 'repeat(' + length + ', minmax(0, 1fr))'
  } else {
    return 'grid-template-columns:' + 'repeat(8, minmax(0, 1fr))'
  }
}

// ---------------------------------- 多实验适配 ----------------------------------

// index 进度条配置
const maxInnerIndex = computed(() => {
  let tempLength = 1
  for (const exp in stepsData[currentIndex.value]) {
    const l = stepsData[currentIndex.value][exp].length
    if (l > maxInnerIndex.value) tempLength = l
  }
  return tempLength
})
const minInnerIndex = ref(1)
const currentInnerIndex = ref(minInnerIndex.value)

// 同时获取多个实验的图像数据
const getMultiImagesData = async (stepData) => {
  const promises = []
  for (const exp in stepData) {
    for (const { filename } of stepData[exp]) {
      if (imagesData[filename]) continue
      // 没有缓存，发起请求
      promises.push(
        new Promise((resolve) => {
          console.log(run_id.value[exp], props.title, filename)
          UTILS.media.get(filename, run_id.value[exp], props.title).then((blob) => {
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
}

const visiableSources = computed(() => {
  const temp = []
  for (const key in stepsData[currentIndex.value]) {
    temp.push(key)
  }
  return temp
})

// ---------------------------------- 数据格式化 ----------------------------------
// 前端映射数据，不包含图像，数据格式：{step: {tag: [{filename: string, caption: string}]}}
const stepsData = {}
// 图像数据缓存, 数据格式：{String<filename>: {blob: Blob, url: string<base64>}}
const imagesData = {}

/**
 * 获取图像数据
 * step 发生改变后触发获取
 */
const getImagesData = async (stepData) => {
  if (!isMulti.value) return await getSingleImageData(stepData)
  await getMultiImagesData(stepData)
}

/**
 * 单数据源
 */
const getSingleImageData = async (stepData) => {
  const tag = source.value[0]
  const promises = []
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
}

/**
 * 防抖函数，防止请求过于频繁
 */
const debounceGetImagesData = debounce(getImagesData, 500)

/**
 * 将后端传回的数据转换为图像数据，存储到imageStepsData中
 * 后端返回的数据格式为：{exp_name: {..., list:[{data: string or list<strinf>, more: object or list<object>, index: string<number>, create_time: string}]}}
 * 主要工作：
 * 1. 将有用数据按 step 提取到 stepsData
 * 2. 找出最大最小 step
 *
 * @param { Object } data 后端返回的数据
 */
const changeData2Image = (data) => {
  // 遍历数据，将数据转换为图像数据
  let _maxIndex = 0
  let _minIndex = Infinity
  for (const exp in data) {
    if (data[exp] === null) continue
    // 遍历实验下的数据
    _maxIndex = Math.max(Number(data[exp].list[data[exp].list.length - 1].index), _maxIndex)
    _minIndex = Math.min(Number(data[exp].list[0].index), _minIndex)
    for (const item of data[exp].list) {
      // 如果不存在当前 step
      if (!stepsData[item.index]) stepsData[item.index] = {}
      // 对当前实验下的数据进行处理,向stepsData中添加数据, 初始化 step 下的实验存储
      if (!stepsData[item.index][exp]) stepsData[item.index][exp] = []
      // 添加数据,如果data是字符串，则直接添加，如果是数组，则遍历添加
      if (typeof item.data === 'string') {
        stepsData[item.index][exp].push({ filename: item.data, caption: item.more?.caption })
      } else {
        for (let i = 0; i < item.data.length; i++) {
          stepsData[item.index][exp].push({ filename: item.data[i], caption: item.more[i]?.caption })
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

// ---------------------------------- 点击某个图像，放大 ----------------------------------
const isSingleZoom = ref(false)
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
    }

    .image-detail {
      @apply flex flex-col items-center justify-center relative;
    }
  }
}

.download-button {
  @apply opacity-75 absolute top-1 right-1 hidden;
  width: 1.25rem !important;
  height: 1.25rem !important;
  padding: 0.125rem !important;
}

.image-content-no-zoom {
  @apply h-56 overflow-y-auto overflow-x-clip;
}

.image-content-zoom {
  @apply h-[calc(100vh-19rem)] overflow-y-auto overflow-x-clip;
}

.image-single-zoom {
  @apply flex items-center;
  &:hover .download-button {
    @apply block;
  }
}
</style>
