<template>
  <!-- 图表标题 -->
  <p class="text-center font-semibold">{{ title }}</p>
  <!-- 如果图表数据错误 -->
  <div class="flex flex-col justify-center grow text-dimmer gap-2" v-if="error">
    <SLIcon class="mx-auto h-5 w-5" icon="error" />
    <p class="text-center text-xs">
      <!-- 在此处显示错误信息 -->
      <span v-if="!isMulti">{{ $t('chart.charts.image.error', { type: error['data_class'], tag: source[0] }) }}</span>
      <span v-else>{{ $t('chart.error') }}</span>
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
        <ImageModule
          v-for="(s, index) in stepsData[currentIndex][source[0]]"
          :key="index"
          :filename="s.filename"
          :imagesData="imagesData[s.filename].url"
          :caption="s.caption"
          :index="index"
          @zoom="handleClickZoom"
          @download="download"
        />
      </div>
      <!-- 多实验图表 -->
      <div v-if="!loading && isMulti" class="images-container" :style="setGrid(visiableSources.length)">
        <!-- 多实验图表 -->
        <ImageModule
          v-for="(s, name) in stepsData[currentIndex]"
          :key="name"
          :filename="s[currentInnerIndex].filename"
          :imagesData="imagesData[s[currentInnerIndex].filename].url"
          :caption="s[currentInnerIndex].caption"
          :index="name"
          :name="name"
          :color="getColor(name, undefined, chart.id)"
          multi
          @zoom="handleClickZoom"
          @download="download"
        />
      </div>
    </div>
    <div class="md:h-8 md:flex md:gap-10 items-center justify-center mt-2">
      <SlideBar
        :class="{ 'md:!justify-end': isMulti }"
        v-model="currentIndex"
        :max="maxIndex"
        :min="minIndex"
        :bar-color="barColor"
        :key="slideKey"
        @turn="handleTurn"
        v-if="maxIndex !== minIndex"
      />
      <SlideBar
        :class="{ 'md:!justify-start': maxIndex !== minIndex }"
        :style="{ visibility: maxInnerIndex || 'hidden' }"
        v-model="currentInnerIndex"
        :max="maxInnerIndex"
        :min="minInnerIndex"
        :bar-color="barColor"
        :key="maxInnerIndex"
        reference="Index"
        @turn="handleTurnIndex"
        v-if="isMulti"
      />
    </div>
    <!-- 放大效果弹窗 -->
    <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom" @onExit="exitByEsc">
      <p class="text-center mt-4 mb-10 text-2xl font-semibold">{{ title }}</p>
      <div class="image-content image-content-zoom">
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
          <ImageModule
            v-for="(s, index) in stepsData[currentIndex][source[0]]"
            :key="index"
            :filename="s.filename"
            :imagesData="imagesData[s.filename].url"
            :caption="s.caption"
            :index="index"
            @zoom="handleClickZoom"
            @download="download"
          />
        </div>
        <!-- 多实验图表 -->
        <div v-if="!loading && isMulti" class="images-container" :style="setGrid(visiableSources.length)">
          <!-- 多实验图表 -->
          <ImageModule
            v-for="(s, name) in stepsData[currentIndex]"
            :key="name"
            :filename="s[currentInnerIndex].filename"
            :imagesData="imagesData[s[currentInnerIndex].filename].url"
            :caption="s[currentInnerIndex].caption"
            :index="name"
            :name="name"
            :color="getColor(name, undefined, chart.id)"
            multi
            @zoom="handleClickZoom"
            @download="download"
          />
        </div>
      </div>
      <div class="md:h-8 md:flex md:gap-10 items-center justify-center mt-2">
        <SlideBar
          :class="{ 'md:!justify-end': isMulti }"
          v-model="currentIndex"
          :max="maxIndex"
          :min="minIndex"
          :bar-color="barColor"
          :key="slideKey"
          @turn="handleTurn"
          v-if="maxIndex !== minIndex"
          :turn-by-arrow="isZoom && !isSingleZoom"
        />
        <SlideBar
          :class="{ 'md:!justify-start': maxIndex !== minIndex }"
          :style="{ visibility: maxInnerIndex || 'hidden' }"
          v-model="currentInnerIndex"
          :max="maxInnerIndex"
          :min="minInnerIndex"
          :bar-color="barColor"
          :key="maxInnerIndex"
          reference="index"
          @turn="handleTurnIndex"
          v-if="isMulti"
          :turn-by-arrow="isZoom && !isSingleZoom && isMulti"
          vertical-arrow
        />
      </div>
    </SLModal>
    <!-- 额外的放大功能，点击某个图像，放大显示 -->
    <SLModal
      class="w-full flex justify-center min-h-[calc(100vh-8rem)] p-14 relative"
      max-w="-1"
      v-model="isSingleZoom"
      close-on-overlay-click
      @onExit="exitByEsc"
    >
      <!-- 标题 -->
      <p class="w-full overflow-hidden text-center text-lg font-semibold absolute top-5">
        <span v-if="isMulti">{{ visiableSources[currentSingleImageIndex] }} / </span>
        {{ signleZoomFilename }}
      </p>
      <!-- 图片 -->
      <div class="image-single-zoom">
        <!-- 上一张图片 -->
        <SLIcon icon="down" class="icon rotate-90" @click="handleSingleChange({ key: 'ArrowLeft' })"></SLIcon>
        <!-- 当前图片 -->
        <div class="relative mx-5 select-none">
          <img :src="imagesData[signleZoomFilename].url" class="w-full" />
          <DownloadButton class="image-download-button" @click.stop="download(signleZoomFilename)" />
        </div>
        <!-- 下一张图片 -->
        <SLIcon icon="down" class="icon -rotate-90" @click="handleSingleChange({ key: 'ArrowRight' })"></SLIcon>
      </div>
      <p class="w-full text-center absolute bottom-5 select-none">
        {{
          `${currentSingleImageIndex + 1} / ${
            isMulti ? visiableSources.length : stepsData[currentIndex][source[0]].length
          }`
        }}
      </p>
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
import SlideBar from '../components/SlideBar.vue'
import { debounce } from '@swanlab-vue/utils/common'
import DownloadButton from '../components/DownloadButton.vue'
import ImageModule from '../modules/ImageModule.vue'

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
const media = inject('media')
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
  return props.chart.multi
})

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------

const error = computed(() => {
  if (!props.chart.error) return false
  if (!isMulti.value) return props.chart.error[source.value[0]]
  // 如果是多实验，检查每个实验的error
  for (const exp of source.value) {
    // 如果有错误，直接返回
    if (props.chart.error[exp]) return true
  }
  return false
})

// ---------------------------------- 图表颜色配置 ----------------------------------
// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const getColor = inject('getColor')
const defaultColor = inject('defaultColor')
// ---------------------------------- 组件渲染逻辑 ----------------------------------
// 已经滑动部分颜色，应该通过色盘计算得到
const barColor = defaultColor
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
    if (val in stepsData.value) {
      if (val === __currentIndex.value && loading.value === false) return
      __currentIndex.value = val
    } else {
      // 寻找一个最接近的值
      const keys = Array.from(Object.keys(stepsData.value))
      const index = keys.findIndex((item) => item > val)
      const num = Number(index === -1 ? keys[keys.length - 1] : keys[index])
      if (num === __currentIndex.value) return
      __currentIndex.value = num
    }
    loading.value = true
    imageContentRef.value.height = imageContentRef.value.offsetHeight + 'px'
    // 如果是多实验模式，将内部数据索引置 0
    if (isMulti.value) currentInnerIndex.value = 0
    // 获取数据
    debounceGetImagesData(stepsData.value[__currentIndex.value])
  }
})

// 事件处理，触发slideBar的turn事件
const handleTurn = (direction, value) => {
  const keys = Array.from(Object.keys(stepsData.value))
  const index = keys.findIndex((item) => item > value)
  if (direction === 'forward') {
    currentIndex.value = Number(index === -1 ? keys[keys.length - 1] : keys[index])
  } else {
    // 向下获取
    if (index === -1) currentIndex.value = Number(keys[Math.max(0, keys.length - 2)])
    else currentIndex.value = Number(keys[Math.max(0, index - 2)])
  }
}

// 在多实验图表时完成 index 的翻页
const handleTurnIndex = (direction, value) => {
  currentInnerIndex.value = direction === 'forward' ? value + 1 : value - 1
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
  let tempLength = 0
  for (const exp in stepsData.value[currentIndex.value]) {
    const l = stepsData.value[currentIndex.value][exp].length - 1
    if (l > tempLength) tempLength = l
  }
  return tempLength
})
const minInnerIndex = ref(0)
const currentInnerIndex = ref(minInnerIndex.value)

// 同时获取多个实验的图像数据
const getMultiImagesData = async (stepData) => {
  const promises = []
  for (const exp in stepData) {
    for (const { filename, experiment_id } of stepData[exp]) {
      if (imagesData[filename]) continue
      // 没有缓存，发起请求
      promises.push(
        new Promise((resolve) => {
          media.get(filename, experiment_id, props.title).then((blob) => resolve(transformBlob(blob, filename)))
        })
      )
    }
  }

  await Promise.all(promises)
}

/**
 * 多实验图表时，当前 step 中含有的实验名数组
 */
const visiableSources = computed(() => {
  const temp = []
  for (const key in stepsData.value[currentIndex.value]) {
    temp.push(key)
  }
  return temp
})

// ---------------------------------- 数据格式化 ----------------------------------
// 前端映射数据，不包含图像，数据格式：{step: {tag: [{filename: string, caption: string}]}}
const stepsData = ref({})
// 图像数据缓存, 数据格式：{String<filename>: {blob: Blob, url: string<base64>}}
const imagesData = {}

/**
 * 获取图像数据
 * step 发生改变后触发获取
 */
const getImagesData = async (stepData) => {
  if (!isMulti.value) await getSingleImageData(stepData)
  else await getMultiImagesData(stepData)
  loading.value = false
  imageContentRef.value.height = ''
}

/**
 * 单数据源
 */
const getSingleImageData = async (stepData) => {
  const tag = source.value[0]
  const promises = []
  for (const { filename, experiment_id } of stepData[tag]) {
    if (!imagesData[filename]) {
      promises.push(
        new Promise((resolve) => {
          media.get(filename, experiment_id, tag).then((blob) => resolve(transformBlob(blob, filename)))
        })
      )
    }
  }
  await Promise.all(promises)
}

/**
 * 将 blob 转成图片
 */
const transformBlob = (blob, filename) => {
  return new Promise((resolve) => {
    // blob转换为图像base64
    const reader = new FileReader()
    reader.onload = function () {
      imagesData[filename] = { blob, url: reader.result }
      resolve()
    }
    reader.readAsDataURL(blob)
  })
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
  for (const source in data) {
    if (data[source] === null) continue
    // 遍历实验下的数据
    _maxIndex = Math.max(Number(data[source].list[data[source].list.length - 1].index), _maxIndex)
    _minIndex = Math.min(Number(data[source].list[0].index), _minIndex)
    const experiment_id = props.chart.source_map[source]
    for (const item of data[source].list) {
      // 如果不存在当前 step
      if (!stepsData.value[item.index]) stepsData.value[item.index] = {}
      // 对当前实验下的数据进行处理,向stepsData中添加数据, 初始化 step 下的实验存储
      if (stepsData.value[item.index][source]) continue
      stepsData.value[item.index][source] = []
      // 添加数据,如果data是字符串，则直接添加，如果是数组，则遍历添加
      if (typeof item.data === 'string') {
        stepsData.value[item.index][source].push({ filename: item.data, caption: item.more?.caption, experiment_id })
      } else {
        for (let i = 0; i < item.data.length; i++) {
          stepsData.value[item.index][source].push({
            filename: item.data[i],
            caption: item.more[i]?.caption,
            experiment_id
          })
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
  currentIndex.value = maxIndex.value
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
// 当前单个图像的索引
const currentSingleImageIndex = ref(0)
// 点击某个图像，放大
const handleClickZoom = (filename, index) => {
  signleZoomFilename.value = filename
  isSingleZoom.value = true
  currentSingleImageIndex.value = isMulti.value ? visiableSources.value.indexOf(index) : index
}

// ---------------------------------- ESC 退出弹窗 ----------------------------------

// 通过 esc 按键关闭弹窗
const exitByEsc = () => {
  // 如果两个弹窗都有，只关闭单图弹窗
  if (isZoom.value && isSingleZoom.value) return (isSingleZoom.value = false)
  // 剩下的情况都是只有一个弹窗，全部关闭
  isZoom.value = false
  isSingleZoom.value = false
}

// ---------------------------------- 左右方向键翻页 ----------------------------------

// 方向键切换单图 - 回调
const handleSingleChange = ({ key }) => {
  const isSingle = !isMulti.value
  const images = isSingle ? stepsData.value[currentIndex.value][source.value[0]] : stepsData.value[currentIndex.value]

  if (
    (key === 'ArrowRight' &&
      currentSingleImageIndex.value < (isSingle ? images.length - 1 : visiableSources.value.length - 1)) ||
    (key === 'ArrowLeft' && currentSingleImageIndex.value > 0)
  ) {
    currentSingleImageIndex.value += key === 'ArrowRight' ? 1 : -1
  }

  const filename = isSingle
    ? images[currentSingleImageIndex.value].filename
    : images[visiableSources.value[currentSingleImageIndex.value]][currentInnerIndex.value].filename

  signleZoomFilename.value = filename
}

// 订阅通过左右方向键切换单图
watch(
  () => isSingleZoom.value,
  (newVal) => {
    if (newVal) {
      window.addEventListener('keyup', handleSingleChange)
    } else {
      window.removeEventListener('keyup', handleSingleChange)
    }
  }
)

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
  .images-container {
    @apply grid gap-3 h-full;
  }
}
.image-content-no-zoom {
  @apply h-56 overflow-y-auto overflow-x-clip;
}

.image-content-zoom {
  @apply h-[calc(100vh-19rem)] overflow-y-auto overflow-x-clip;
}

.image-single-zoom {
  @apply flex items-center justify-between w-full;
  &:hover .image-download-button {
    @apply block;
  }

  .icon {
    @apply w-10 h-10 cursor-pointer border rounded-full opacity-20 transition-all;

    &:hover {
      @apply opacity-100;
    }
  }
}
</style>

<style lang="scss">
.image-download-button {
  @apply opacity-75 absolute top-1 right-1 hidden;
  width: 1.25rem !important;
  height: 1.25rem !important;
  padding: 0.125rem !important;
}
</style>
