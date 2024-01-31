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
      <AudioModule :audios="audioData" :key="nowStep" v-if="audioData" />
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
const minIndex = ref(undefined)
// 最大索引
const maxIndex = ref(undefined)
// slide的key
const slideKey = ref(0)
watch([maxIndex, minIndex], ([max, min]) => {
  slideKey.value = max + '-' + min
})

// 已经滑动部分颜色，应该通过色盘计算得到
const barColor = inject('colors')[0]

// ---------------------------------- 数据驱动控制音频组件渲染 ----------------------------------
const audiosTagData = ref()
// 音频数据缓存，格式为 {step: nowData }
const audiosData = ref({})
// 当前展示的是第几步的数据，格式为{tag: {audio: AudioBuffer, title: String}}
const audioData = computed(() => {
  if (nowStep.value === undefined) return undefined
  // console.log('nowStep', nowStep.value, audiosData.value[nowStep.value])
  return audiosData.value[nowStep.value]
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

/**
 * 获取音频数据
 * @param {str} filename 文件名
 */
const getMediaData = async (data, tag) => {
  const audioBlob = await UTILS.media.get(data, run_id.value, tag)
  // console.log('audioBlob', audioBlob)
  // 将blob转换成AudioBuffer
  // 1. 将 Blob 转换为 ArrayBuffer
  const arrayBuffer = await audioBlob.arrayBuffer()

  // 2. 使用 Web Audio API 的 decodeAudioData 方法将 ArrayBuffer 转换为 AudioBuffer
  const audioContext = new (window.AudioContext || window.webkitAudioContext)()

  return new Promise((resolve, reject) => {
    audioContext.decodeAudioData(
      arrayBuffer,
      (buffer) => {
        resolve(buffer)
      },
      (error) => {
        reject(error)
      }
    )
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

/**
 * 将 tag 数据转成 AudioBuffer，添加到audiosData中
 * 并且设置当前step为index对应的step,最终返回当前的step和最大step和最小step
 * 后续如果要兼容多数据，就在此处添加逻辑
 * @param { number } index
 */
const tagData2Buffer = async (index) => {
  const tag = source.value[0]
  const data = audiosTagData.value[tag].list[index]
  const audioBuffer = await getMediaData(data.data, tag)
  nowStep.value = Number(data.index)
  return {
    nowStep: Number(data.index),
    maxStep: Number(audiosTagData.value[tag].list[audiosTagData.value[tag].list.length - 1].index),
    minStep: Number(audiosTagData.value[tag].list[0].index),
    data: [{ audioBuffer, title: data.data, tag: tag }],
    tag
  }
}

// ---------------------------------- 渲染、重渲染功能 ----------------------------------
// 渲染
const render = async (data) => {
  audiosTagData.value = data
  // 根据 data 获取所有的音频数据，存储到audiosData中
  // 当前step必然为第一个，因为是第一次渲染，所以其他的数据不需要请求,只需要请求第一个就可以
  const steps = await tagData2Buffer(0)
  // 将数据添加到audiosData中
  audiosData.value[steps.nowStep] = steps.data
  // console.log('steps', steps, audiosData.value)
  // 设置最大值最小值
  currentIndex.value = nowStep.value = steps.nowStep
  maxIndex.value = steps.maxStep
  minIndex.value = steps.minStep
}

// 重渲染，只需要往audiosData中添加数据并更新maxIndex和minIndex即可，不需要重复请求数据
const change = async (data) => {
  audiosTagData.value = data
  // 设置最大值最小值
  const tag = source.value[0]
  maxIndex.value = Number(data[tag].list[data[tag].list.length - 1].index)
  minIndex.value = Number(data[tag].list[0].index)
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
