<template>
  <!-- 图表标题 -->
  <p class="text-center font-semibold">{{ title }}</p>
  <SlideBar v-model="now" :max="max" :min="min" :bar-color="barColor" />
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
    <div ref="canvas"></div>
    <!-- 放大效果弹窗 -->
    <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom"> </SLModal>
  </template>
  {{ display }}
</template>

<script setup>
/**
 * @description:
 * @file: AudioChart.vue
 * @since: 2024-01-29 20:31:18
 **/
import SLModal from '@swanlab-vue/components/SLModal.vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import { onMounted, ref, inject, computed } from 'vue'
import SlideBar from '../components/SlideBar.vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import * as UTILS from './utils'
import http from '@swanlab-vue/api/http'
import { useExperimentStore } from '@swanlab-vue/store'

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
  return props.chart?.source
})

// 音频上下文
let audioCtx = new (window.AudioContext || window.webkitAudioContext)()
// let audioSource = audioCtx.createBufferSource()

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------

const error = ref(props.chart.error)

// ---------------------------------- 图表颜色配置 ----------------------------------

// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const colors = inject('colors')
if (!colors) throw new Error('colors is not defined, please provide colors in parent component')

// ---------------------------------- 组件渲染逻辑 ----------------------------------

// 当前展示的是对应 tag 的第几条数据
const display = ref({})

// 初始化显示配置
onMounted(() => {
  // 每一个相关 tag 都展示第一条数据
  source.value.forEach((tag) => {
    display.value[tag] = 1
  })
})

// 画布元素
const canvas = ref(null)

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
 * @param {*} tag
 */
const format = async (tag) => {
  let promises = []
  // 将 blob 转成 ArrayBuffer，因为 decodeAudioData 接收参数类型为 ArrayBuffer，因为
  audioBlob.value[tag].forEach((blob) => {
    promises.push(blobToBuf(blob))
  })
  /**
   * 这里需要注意：
   * buffers 中的元素为对象
   * - buf : ArrayBuffer
   * - type: BlobType
   */
  const buffers = await Promise.all(promises)
  promises = []
  // 将所有缓存都转成 AudioBuffer
  buffers.forEach((buffer) => {
    promises.push(audioCtx.decodeAudioData(buffer.buf))
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
  console.log(audioData.value)
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
        var buf = e.target.result
        resolve({ buf, type })
      },
      false
    )
  })
}

// ---------------------------------- 渲染、重渲染功能 ----------------------------------

// 渲染
const render = async (data) => {
  if (source.value && source.value.length < 0) {
    error.value = true
    return
  }
  await getMediaData(data)
  Object.keys(audioBlob.value).forEach(format)
}
// 重渲染
const change = (data) => {}

// ---------------------------------- 放大功能 ----------------------------------
// // 是否放大
// const isZoom = ref(false)
// // 放大数据
// const zoom = (data) => {
//   isZoom.value = true
//   // 放大后图表的高度
//   const height = window.innerHeight * 0.6
//   addTaskToBrowserMainThread(() => {})
// }

// ---------------------------------- 实例：滑块的使用 ----------------------------------
// 当前值
const now = ref(0)
// 最大值
const max = ref(100)
// 最小值
const min = ref(0)
// 已经滑动部分颜色，应该通过色盘计算得到
const barColor = inject('colors')[0]

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  render,
  change
  // zoom
})
</script>

<style lang="scss" scoped></style>
