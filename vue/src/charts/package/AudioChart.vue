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
    <div class="audio-content audio-content-no-zomm" ref="audioContentRef">
      <AudioModule :audios="audioData" :key="currentIndex" v-if="audioData && !loading" />
      <div class="flex flex-col justify-center items-center h-full" v-if="loading">
        <SLLoading />
      </div>
    </div>
    <!-- 滑块 -->
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
        @turn="handleTurnIndex"
        v-if="isMulti"
      />
    </div>
    <!-- 放大效果弹窗 -->
    <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom"> </SLModal>
  </template>
</template>

<script setup>
/**
 * @description:
 * @file: AudioChart.vue
 * @since: 2024-03-09 19:22:22
 **/
import SLModal from '@swanlab-vue/components/SLModal.vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import SlideBar from '../components/SlideBar.vue'
import AudioModule from '../modules/AudioModule.vue'
import { ref, inject } from 'vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import * as UTILS from './utils'
import { useExperimentStore, useProjectStore } from '@swanlab-vue/store'
import { useRoute } from 'vue-router'
import { debounce } from '@swanlab-vue/utils/common'

// ---------------------------------- 配置 ----------------------------------

const experimentStore = useExperimentStore()
const projectStore = useProjectStore()
const route = useRoute()

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

const audioContentRef = ref(null)

/**
 * 图表数据来源，为数组
 * 当isMulti为true时，数组元素为exp，即给该表提供数据的实验的名称
 * 当isMulti为false时，数组元素为tag，即给该表提供数据的tag
 */
const sources = computed(() => {
  return props.chart?.source || []
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
  for (const exp of sources.value) {
    run_ids[exp] = projectStore.getExpRunIdByName(exp)
  }
  return run_ids
})

/**
 * 展示用的音频数据
 */
const audioData = computed(() => {
  const stepData = stepsData[currentIndex.value]
  if (!stepData) return []
  // 如果是单源
  if (!isMulti.value) {
    return stepData[sources.value[0]].map(({ filename, caption }) => {
      return {
        audioBuffer: audiosData[filename],
        title: filename,
        caption
      }
    })
  }
  // 如果是多源
  const data = []
  for (const exp in stepData) {
    const temp = stepData[exp][currentInnerIndex.value]
    data.push({
      audioBuffer: audiosData[temp.filename],
      title: temp.filename,
      caption: temp.caption,
      color: projectStore.colorMap[exp],
      exp
    })
  }
  return data
})

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------

const error = ref(props.chart.error[sources.value[0]])

// ---------------------------------- 图表颜色配置 ----------------------------------

// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const colors = inject('colors')
if (!colors) throw new Error('colors is not defined, please provide colors in parent component')

// ---------------------------------- step 滑块 ----------------------------------

// 已经滑动部分颜色，应该通过色盘计算得到
const barColor = inject('colors')[0]
// 当前滑块索引
const __currentIndex = ref(0)
// 最小索引
const minIndex = ref(NaN)
// 最大索引
const maxIndex = ref(NaN)
// slide的key
const slideKey = ref(0)
watch([maxIndex, minIndex], ([max, min]) => {
  slideKey.value = max + '-' + min
})
const loading = ref(true)

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
    // audioContentRef.value.style.height = audioContentRef.value.offsetHeight + 'px'
    loading.value = true
    // 如果是多实验模式，将内部数据索引置 0
    if (isMulti.value) currentInnerIndex.value = 0
    // 获取数据
    debounceGetAudiosData(stepsData[__currentIndex.value])
  }
})

/**
 * step 滑块，点击上下按钮翻页
 */
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

// ---------------------------------- 数据格式化 ----------------------------------

// 关联 step 与数据
const stepsData = {}
// 缓存
const audiosData = reactive({})

/**
 * 解析 log 数据
 */
const changeData2Audio = (data) => {
  let _maxIndex = 0
  let _minIndex = Infinity
  for (const source in data) {
    if (data[source] === null) continue
    // 获取最大和最小索引
    _maxIndex = Math.max(Number(data[source].list[data[source].list.length - 1].index), _maxIndex)
    _minIndex = Math.min(Number(data[source].list[0].index), _minIndex)

    for (const item of data[source].list) {
      // 如果不存在当前 step
      if (!stepsData[item.index]) stepsData[item.index] = {}
      // 对当前实验下的数据进行处理,向stepsData中添加数据, 初始化 step 下的实验存储
      if (!stepsData[item.index][source]) stepsData[item.index][source] = []
      // 添加数据,如果data是字符串，则直接添加，如果是数组，则遍历添加
      if (typeof item.data === 'string') {
        stepsData[item.index][source].push({ filename: item.data, caption: item.more?.caption })
      } else {
        for (let i = 0; i < item.data.length; i++) {
          stepsData[item.index][source].push({ filename: item.data[i], caption: item.more[i]?.caption })
        }
      }
    }
  }
  maxIndex.value = _maxIndex
  minIndex.value = _minIndex
}

// ---------------------------------- 数据请求和解析 ----------------------------------

const debounceGetAudiosData = debounce(async (stepData) => {
  loading.value = true

  const promises = []
  // 遍历该 step 下的数据源
  for (const source in stepData) {
    // 如果数据源中没有东西，跳过
    if (stepData[source] === null) continue
    // 遍历该数据源下的数据
    for (const { filename } of stepData[source]) {
      if (audiosData[filename]) continue
      // 没有缓存，需要请求
      promises.push(
        new Promise((resolve) => {
          const runid = isMulti.value ? run_id.value[source] : run_id.value
          UTILS.media.get(filename, runid, props.title).then((blob) => resolve(transformBlob(blob, filename)))
        })
      )
    }
  }

  await Promise.all(promises)
  loading.value = false
}, 500)

/**
 * 将 blob 转成 AudioBuffer
 */
const transformBlob = async (blob, filename) => {
  const arrayBuffer = await blob.arrayBuffer()
  const audioContext = new (window.AudioContext || window.webkitAudioContext)()
  return new Promise((resolve, reject) => {
    audioContext.decodeAudioData(
      arrayBuffer,
      (audioBuffer) => {
        audiosData[filename] = audioBuffer
        resolve(audioBuffer)
      },
      (error) => {
        error.value = error
        reject(error)
      }
    )
  })
}

// ---------------------------------- 多实验 ----------------------------------

// index 进度条配置
const maxInnerIndex = computed(() => {
  let tempLength = 0
  for (const exp in stepsData[currentIndex.value]) {
    const l = stepsData[currentIndex.value][exp].length - 1
    if (l >= maxInnerIndex.value) tempLength = l
  }
  return tempLength
})
const minInnerIndex = ref(0)
const currentInnerIndex = ref(minInnerIndex.value)

const handleTurnIndex = (direction, value) => {
  currentInnerIndex.value = direction === 'forward' ? value + 1 : value - 1
}

// ---------------------------------- 渲染、重渲染功能 ----------------------------------

// 渲染
const render = (data) => {
  changeData2Audio(data)
  currentIndex.value = minIndex.value
}
// 重渲染
const change = (data) => {
  changeData2Audio(data)
}

// ---------------------------------- 放大功能 ----------------------------------
// 是否放大
const isZoom = ref(false)
// 放大数据
const zoom = (data) => {
  isZoom.value = true
  // 放大后图表的高度
  const height = window.innerHeight * 0.6
  addTaskToBrowserMainThread(() => {})
}

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  render,
  change,
  zoom
})
</script>

<style lang="scss" scoped></style>
