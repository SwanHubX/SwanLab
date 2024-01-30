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
    <!-- 放大效果弹窗 -->
    <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom"> </SLModal>
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
import { ref, inject, computed } from 'vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import * as UTILS from './utils'
import http from '@swanlab-vue/api/http'
import { useExperimentStore } from '@swanlab-vue/store'

const experimentStore = useExperimentStore()
const run_id = computed(() => experimentStore.experiment.run_id)

let audioCtx = new (window.AudioContext || window.webkitAudioContext)()
// let audioSource = audioCtx.createBufferSource()
// ---------------------------------- 配置 ----------------------------------
const props = defineProps({
  title: {
    type: String,
    required: true
  },
  chart: {
    type: Object,
    required: true
  }
})

// 图表相关 tag
const source = computed(() => {
  return props.chart?.source
})

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------

const error = ref(props.chart.error)

// ---------------------------------- 图表颜色配置 ----------------------------------
// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const colors = inject('colors')
if (!colors) throw new Error('colors is not defined, please provide colors in parent component')
// ---------------------------------- 组件渲染逻辑 ----------------------------------

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
  res.forEach((item, index) => {
    audioBlob.value[tags[index]] = item
  })
}

const getMedia = (tag, list) => {
  const promises = []
  list.forEach((item) => {
    console.log(item)
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

const format = (key) => {
  audioBlob.value[key].forEach((blob) => {
    const list = []
    blobToBuf(blob, (buf, type) => {
      console.log(type)
      audioCtx.decodeAudioData(buf).then((res) => {
        list.push(res)
      })
    })
    console.log(list)
  })
}

const blobToBuf = (blob, callback) => {
  var fr = new FileReader()
  var type = blob.type

  fr.readAsArrayBuffer(blob)
  fr.addEventListener(
    'loadend',
    (e) => {
      var buf = e.target.result
      callback(buf, type)
    },
    false
  )
}

// ---------------------------------- 渲染、重渲染功能 ----------------------------------

// 渲染
const render = async (data) => {
  if (source.value && source.value.length < 0) {
    error.value = true
    return
  }
  await getMediaData(data)
  console.log(audioBlob.value)
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

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  render,
  change
  // zoom
})
</script>

<style lang="scss" scoped></style>
