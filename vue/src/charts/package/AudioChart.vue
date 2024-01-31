<template>
  <!-- 图表标题 -->
  <p class="text-center font-semibold">{{ title }}</p>
  <!-- 如果图表数据错误 -->
  <div class="flex flex-col justify-center grow text-dimmer gap-2" v-if="error">
    <SLIcon class="mx-auto h-5 w-5" icon="error" />
    <p class="text-center text-xs">
      <!-- 在此处显示错误信息 -->
    </p>
  </div>
  <!-- 如果图表数据正确 -->
  <template v-else>
    <!-- 在此处完成图表主体定义 -->
    <div class="mt-1 p-2 w-full border border-dimmer rounded-sm relative h-56">
      <AudioModule :audios="nowData" v-if="nowData" :key="nowStep" />
    </div>
    <SlideBar
      class="mt-2"
      v-model="currentIndex"
      :max="maxIndex"
      :min="minIndex"
      :bar-color="barColor"
      :key="maxIndex"
    />
    <!-- 放大效果弹窗 -->
    <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom">
      <div ref="modalAudioRef"></div>
      <PlayButton @click="playAudio" />
    </SLModal>
  </template>
</template>

<script setup>
/**
 * @description:
 * @file: AudioChart.vue
 * @since: 2024-01-29 20:31:18
 **/
import SLModal from '@swanlab-vue/components/SLModal.vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import { onMounted, watch, ref, inject, computed } from 'vue'
import SlideBar from '../components/SlideBar.vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import * as UTILS from './utils'
import { useExperimentStore } from '@swanlab-vue/store'
import PlayButton from '../components/PlayButton.vue'
import AudioModule from '../modules/AudioModule.vue'

// ---------------------------------- 配置 ----------------------------------

const experimentStore = useExperimentStore()
const run_id = computed(() => experimentStore.experiment.run_id)

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
    default: 0
  }
})

// 图表相关 tag
const source = computed(() => {
  return props.chart?.source || []
})

// 注意！目前是单 tag，故默认选中第一个 tag
const defaultTag = computed(() => {
  return source.value[0]
})

// 音频上下文
const audioCtx = new (window.AudioContext || window.webkitAudioContext)()

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------

const error = ref(props.chart.error)

// ---------------------------------- 图表颜色配置 ----------------------------------

// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const colors = inject('colors')
if (!colors) throw new Error('colors is not defined, please provide colors in parent component')

// ---------------------------------- 实例：滑块的使用 ----------------------------------

// 当前滑块索引
const currentIndex = ref(0)
// 最小索引
const minIndex = 0
// 最大索引
const maxIndex = computed(() => {
  const tag = defaultTag.value
  return audioData.value[tag]?.length - 1
})
// 已经滑动部分颜色，应该通过色盘计算得到
const barColor = inject('colors')[0]

// ---------------------------------- 数据驱动控制音频组件渲染 ----------------------------------
// 音频数据缓存，格式为 {step: nowData }
const data = ref({})
// 当前展示的是第几步的数据，格式为{tag: {audio: AudioBuffer, title: String}}
const nowData = computed(() => {
  if (nowStep.value === undefined) return undefined
  return data.value[nowStep.value]
})
// 当前展示的步数
const nowStep = ref(undefined)
// ---------------------------------- 组件渲染逻辑 ----------------------------------

// 当前展示的是对应 tag 的第几条数据
const display = ref({})

// 初始化显示配置，
onMounted(() => {
  // 当前为单 tag，不需要这里
  // 每一个相关 tag 都展示第一条数据
  source.value.forEach((tag) => {
    display.value[tag] = 1
  })
})

watch(currentIndex, (newV) => {
  draw(defaultTag.value, newV)
})

/**
 * 播放音频
 */
let audioSource = null
const playAudio = () => {
  if (audioSource !== null) {
    audioSource.stop()
  }
  audioSource = audioCtx.createBufferSource()
  audioSource.buffer = audioData.value[defaultTag.value][currentIndex.value].data
  audioSource.connect(audioCtx.destination)
  audioSource.start()
}

const canvasRef = ref(null)

const draw = (tag, index, ref = canvasRef.value[0], height = 200) => {
  // 创建canvas上下文
  const canvas = ref
  const buffer = getAudioBuffer(tag, index)
  const [positives, negatives] = sample(buffer)

  if (canvas.getContext) {
    const width = canvas.offsetWidth
    const ratio = window.devicePixelRatio || 1
    let ctx = canvas.getContext('2d')
    canvas.width = Math.round(width * ratio)
    canvas.height = Math.round(height * ratio)
    canvas.style.width = canvas.width + 'px'
    canvas.style.height = canvas.height + 'px'
    let x = 0
    let y = 100
    // let offset = 0
    ctx.fillStyle = colors[0]
    ctx.beginPath()
    ctx.moveTo(x, y)
    // canvas高度200，横坐标在canvas中点100px的位置，横坐标上方绘制正数据，下方绘制负数据
    // 先从左往右绘制正数据
    // x + 0.5是为了解决canvas 1像素线条模糊的问题
    for (let k = 0; k < width; k++) {
      ctx.lineTo(x + k + 0.5, y - 100 * positives[k])
    }

    // 再从右往左绘制负数据
    for (let l = width - 1; l >= 0; l--) {
      ctx.lineTo(x + l + 0.5, y + 100 * Math.abs(negatives[l]))
    }
    // 填充图形
    ctx.fill()
  }
}

/**
 * 抽样
 * @param {AudioBuffer} buffer 音频缓存
 * @returns {array} [positives, negatives] 正数据和负数据
 */
const sample = (buffer) => {
  let data = []
  let originData = buffer.getChannelData(0)
  // 存储所有的正数据
  let positives = []
  // 存储所有的负数据
  let negatives = []
  const step = 10
  // 先每隔100条数据取1条
  for (let i = 0; i < originData.length; i += step) {
    data.push(originData[i])
  }
  // 再从data中每10条取一个最大值一个最小值
  const range = 5
  for (let j = 0, len = parseInt(data.length / range); j < len; j++) {
    let temp = data.slice(j * 10, (j + 1) * 10)
    positives.push(Math.max.apply(null, temp))
    negatives.push(Math.min.apply(null, temp))
  }
  return [positives, negatives]
}

/**
 * 根据 tag 和索引获取 AudioBuffer
 * @param {string} tag
 * @param {int} index
 */
const getAudioBuffer = (tag, index) => {
  return audioData.value[tag][index].data
}

// ---------------------------------- 数据格式化 ----------------------------------

const audioBlob = ref({})
const audioData = ref({})

/**
 * 获取音频数据
 * @param {object} data 与该表相关的 tag 名为 key,对应的值为数组
 */
const getMediaData = async (data) => {
  const tags = []
  const promises = []
  Object.keys(data).forEach((tag) => {
    tags.push(tag)
    promises.push(
      getMedia(
        tag,
        data[tag].list.map((item) => item.data)
      )
    )
  })
  const res = await Promise.all(promises)
  // 保存所有 blob 数据
  res.forEach((item, index) => {
    audioBlob.value[tags[index]] = item
  })
}

/**
 * 获取某个 tag 下所有的音频数据
 * 返回 Promise.all 对象，在 getMediaData 和别的 tag 一起异步获取
 * @param {string} tag tag 名
 * @param {array} list tag 所属的数据列表
 */
const getMedia = (tag, list) => {
  const promises = []
  list.forEach((item) => {
    promises.push(
      new Promise((resolve) => {
        UTILS.media.get(item, run_id.value, tag).then((res) => {
          resolve(res)
        })
      })
    )
  })
  return Promise.all(promises)
}

/**
 * 数据类型转化
 * 从后端获取的数据皆为 Blob，存放在 audioBlob
 * 需要将二进制音频转成音频上下文可以解析和使用的格式 AudioBuffer
 * 使用 decodeAudioData 可以将 ArrayBuffer 转成 AudioBuffer
 *
 * 注意：
 * 一次 format 只是完成了对一个 tag 数据的转化
 *
 * @param {string} tag log tag name
 */
const format = async (tag, notDraw = false) => {
  let promises = []
  // 将 blob 转成 ArrayBuffer，因为 decodeAudioData 接收参数类型为 ArrayBuffer，因为
  audioBlob.value[tag].forEach((blob) => {
    promises.push(blobToBuf(blob))
  })
  /**
   * 这里需要注意：
   * buffers 中的元素为对象
   * - data : ArrayBuffer
   * - type: BlobType
   */
  const buffers = await Promise.all(promises)
  promises = []
  // 将所有缓存都转成 AudioBuffer
  buffers.forEach((buffer) => {
    promises.push(audioCtx.decodeAudioData(buffer.data))
  })
  // 多线程异步转化
  const audioBuffers = await Promise.all(promises)
  // 保存一下
  audioData.value[tag] = audioBuffers.map((audioBuffer, index) => {
    return {
      data: audioBuffer,
      type: buffers[index].type
    }
  })
  // 到这里所有数据都被转成 AudioBuffer 类型，可以进行抽样和绘制
  // 在这里只绘制当前 tag 的部分
  // 当前只考虑单tag音频图表，用 currentIndex ，若图表需要含有多个 tag 则用 display.value[tag]
  if (notDraw) return // 如果只需要完成格式化，而不需要画表，直接返回
  draw(tag, currentIndex.value)
}

/**
 * 将 blob 数据转成 ArrayBuffer
 * @param {blob} blob
 */
const blobToBuf = (blob) => {
  return new Promise((resolve) => {
    var fr = new FileReader()
    var type = blob.type

    fr.readAsArrayBuffer(blob)
    fr.addEventListener(
      'loadend',
      (e) => {
        var data = e.target.result
        resolve({ data, type })
      },
      false
    )
  })
}

// ---------------------------------- 渲染、重渲染功能 ----------------------------------

const chartData = ref([])

// 渲染
const render = async (data) => {
  chartData.value = data
  // console.log(data)
  if (source.value && source.value.length < 0) {
    error.value = true
    return
  }
  await getMediaData(data)
  // 对全部数据进行格式化
  Object.keys(audioBlob.value).forEach(format)
}
// 重渲染
const change = async (data) => {
  await getMediaData(data)
  // 更新的时候，只需要对 data 中包含的 tag，也就是有改动的 tag 进行格式化
  Object.keys(data).forEach((key) => format(key, true))
}

// ---------------------------------- 放大功能 ----------------------------------
// 是否放大
const isZoom = ref(false)
// 弹窗画板
const modalCanvasRef = ref(null)
// 放大数据
const zoom = () => {
  isZoom.value = true
  addTaskToBrowserMainThread(() => {
    draw(defaultTag.value, currentIndex.value, modalCanvasRef.value)
  })
}

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  render,
  change,
  zoom
})
</script>

<style lang="scss" scoped>
canvas {
  @apply relative;
  &::before {
    @apply absolute left-0 top-0 w-full h-full bg-positive-highest border-x z-10;
    content: '111';
  }
}
</style>
