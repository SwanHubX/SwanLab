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
      <div class="flex flex-col justify-center items-center h-full" v-if="loading">
        <SLLoading />
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
      <div class="mt-15 p-2 w-full border border-dimmer rounded-sm relative h-56">
        <AudioModule :audios="audioData" :key="nowStep" v-if="audioData" />
      </div>
      <div class="h-8 mt-20">
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
import { watch, ref, inject, computed } from 'vue'
import SlideBar from '../components/SlideBar.vue'
import { debounce } from '@swanlab-vue/utils/common'
import * as UTILS from './utils'
import { useExperimentStore } from '@swanlab-vue/store'
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
const loading = ref(false)

// 已经滑动部分颜色，应该通过色盘计算得到
const barColor = inject('colors')[0]

// ---------------------------------- 控制渲染第几个数据 ----------------------------------

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

// ---------------------------------- 数据格式化 ----------------------------------

/**
 * 获取音频数据
 * @param {str} filename 文件名
 * @param {str} tag log tag 名称
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
// 放大数据
const zoom = () => {
  isZoom.value = true
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
