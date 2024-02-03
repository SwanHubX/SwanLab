<template>
  <!-- 图表标题 -->
  <p class="text-center font-semibold">{{ title }}</p>
  <!-- 如果图表数据错误 -->
  <div class="flex flex-col justify-center grow text-dimmer gap-2" v-if="error">
    <SLIcon class="mx-auto h-5 w-5" icon="error" />
    <p class="text-center text-xs">
      <!-- 在此处显示错误信息 -->
      {{ $t('common.chart.charts.audio.error', { type: error['data_class'], tag: source[0] }) }}
    </p>
  </div>
  <!-- 如果图表数据正确 -->
  <template v-else>
    <!-- 在此处完成图表主体定义 -->
    <div class="audio-content" ref="audioContentRef">
      <AudioModule :audios="audioData" :key="nowStep" v-if="audioData && !loading" />
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
      <p class="text-center mt-4 mb-10 text-2xl font-semibold">{{ title }}</p>
      <div class="audio-content" ref="audioContentRef">
        <AudioModule :audios="audioData" :key="nowStep" v-if="audioData && !loading" />
        <div class="flex flex-col justify-center items-center h-full" v-if="loading">
          <SLLoading />
        </div>
      </div>
      <div class="h-8 mt-10">
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

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------

const error = ref(props.chart.error)

// ---------------------------------- 图表颜色配置 ----------------------------------

// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const colors = inject('colors')
if (!colors) throw new Error('colors is not defined, please provide colors in parent component')

// ---------------------------------- 实例：滑块的使用 ----------------------------------

const audioContentRef = ref(null)
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
const loading = ref(false)

const currentIndex = computed({
  get: () => __currentIndex.value,
  set: (val) => {
    if (stepMap.has(val)) {
      // 如果存在，则设置当前值
      if (val === __currentIndex.value) return
      __currentIndex.value = val
      // 当前值
      // console.log('设置当前值为', val)
    } else {
      // console.log('当前值不存在: ', val)
      // 寻找最近的值
      const keys = Array.from(stepMap.keys())
      const index = keys.findIndex((item) => item > val)
      const num = Number(index === -1 ? keys[keys.length - 1] : keys[index])
      if (num === __currentIndex.value) return
      __currentIndex.value = num
    }
    // 设置audioContentRef的height为当前height
    // console.log('设置audioContentRef的height为当前height', audioContentRef.value)
    audioContentRef.value.style.height = audioContentRef.value.offsetHeight + 'px'
    loading.value = true
    debounceRender()
  }
})

// ---------------------------------- 控制渲染第几个数据，防抖 ----------------------------------

const debounceRender = debounce(() => {
  const step = __currentIndex.value

  if (!audiosData.value[step]) {
    tagData2Buffer(stepMap.get(step)[0][1]).then((data) => {
      loading.value = false
      console.log('请求数据成功，设置当前step为', step)
      nowStep.value = step
      // 将数据添加到audiosData中
      audiosData.value[data.nowStep] = data.data
      // 设置audioContentRef的height
      audioContentRef.value.style.height = ''
    })
  } else {
    console.log('已经存在数据，不需要请求, 设置当前step为', step)
    loading.value = false
    nowStep.value = step
    // 设置audioContentRef的height
    audioContentRef.value.style.height = ''
  }
}, 500)

// ---------------------------------- 数据驱动控制音频组件渲染 ----------------------------------
// 音频数据，格式为 {tag: {audio: AudioBuffer, title: String}}
const audiosTagData = ref()
// 音频数据缓存，格式为 {step: nowData }
const audiosData = ref({})
// 当前展示的是第几步的数据，格式为{tag: {audio: AudioBuffer, title: String}}
const audioData = computed(() => {
  if (nowStep.value === undefined) return undefined
  // console.log('nowStep变化：', nowStep.value)
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
      (audioBuffer) => {
        resolve({ audioBuffer, audioBlob })
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
 * @param { number } index 当前数据的索引,并不是step
 */
const tagData2Buffer = async (index) => {
  const tag = source.value[0]
  // console.log('tagData2Buffer', tag, index)
  const data = audiosTagData.value[tag].list[index]
  // data.data可能是list<string>或者string
  const returnData = []
  if (Array.isArray(data.data)) {
    for (let i = 0; i < data.data.length; i++) {
      const item = data.data[i]
      const { audioBuffer, audioBlob } = await getMediaData(item, tag)
      audiosData.value[index] = audiosData.value[index] || []
      returnData.push({ audioBuffer, audioBlob, title: item, tag: tag, caption: data.more[i]?.caption })
    }
  } else {
    const { audioBuffer, audioBlob } = await getMediaData(data.data, tag)
    audiosData.value[index] = audiosData.value[index] || []
    returnData.push({ audioBuffer, audioBlob, title: data.data, tag: tag, caption: data.more?.caption })
  }
  nowStep.value = Number(data.index)
  return {
    nowStep: Number(data.index),
    maxStep: Number(audiosTagData.value[tag].list[audiosTagData.value[tag].list.length - 1].index),
    minStep: Number(audiosTagData.value[tag].list[0].index),
    data: returnData,
    tag
  }
}

//  step与data[tag].list[index]的映射关系，每次change都会更新
// 通过step，可以找到对应的tag和index，{step: [[tag, index],]}
const stepMap = new Map()

const setStepMap = (data) => {
  // 遍历data[tag].list，将step与index的映射关系存储到stepMap中
  // 获取所有tag
  const tags = Object.keys(data)
  // 最长数组长度
  const maxLength = Math.max(...tags.map((tag) => data[tag].list.length))
  for (let i = 0; i < maxLength; i++) {
    for (const tag of tags) {
      const step = Number(data[tag].list[i].index)
      if (!stepMap.has(i)) stepMap.set(step, [])
      if (data[tag].list[i]) {
        // 如果这个step中已经存在此tag，则不需要添加
        if (stepMap.get(step).findIndex((item) => item[0] === tag) === -1) {
          stepMap.get(step).push([tag, i])
        }
      }
    }
  }
}

// ---------------------------------- 渲染、重渲染功能 ----------------------------------
// 渲染
const render = async (data) => {
  audiosTagData.value = data
  setStepMap(data)
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
  setStepMap(data)
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
.audio-content {
  @apply mt-1 p-2 w-full border border-dimmer rounded-sm relative min-h-[224px];
}
</style>
