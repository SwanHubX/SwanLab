<template>
  <!-- 图表标题 -->
  <p class="text-center font-semibold">{{ title }}</p>
  <!-- 如果图表数据错误 -->
  <div class="flex flex-col justify-center grow text-dimmer gap-2" v-if="error">
    <SLIcon class="mx-auto h-5 w-5" icon="error" />
    <p class="text-center text-xs">
      <!-- 在此处显示错误信息 -->
      {{ error }}
    </p>
  </div>
  <!-- 如果图表数据正确 -->
  <template v-else>
    <!-- 在此处完成图表主体定义 -->
    <TextModule :data="original_data" :source="source"></TextModule>
    <!-- 放大效果弹窗 -->
    <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom">
      <TextModule :data="original_data" :source="source" />
    </SLModal>
  </template>
</template>

<script setup>
/**
 * @description: 文字图标
 * @file: TextChart.vue
 * @since: 2024-02-20 20:04:09
 **/
import SLModal from '@swanlab-vue/components/SLModal.vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import { ref, inject, computed } from 'vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import * as UTILS from './utils'
import TextModule from '../modules/TextModule.vue'

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

const original_data = ref()
const source = ref(props.chart.source)

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------

const error = computed(() => {
  if (!props.chart.error) {
    return props.chart.error
  } else if (!original_data.value) {
    return 'No Data'
  }
  return false
})

// ---------------------------------- 图表颜色配置 ----------------------------------
// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const colors = inject('colors')
if (!colors) throw new Error('colors is not defined, please provide colors in parent component')
// ---------------------------------- 组件渲染逻辑 ----------------------------------

// ---------------------------------- 数据格式化 ----------------------------------

// ---------------------------------- 渲染、重渲染功能 ----------------------------------

// 渲染
const render = (data) => {
  console.log('render')
  original_data.value = data
}

// 重渲染
const change = (data) => {
  // 将发生更新的 tag 数据保存到原始数据中
  for (let key in data) {
    original_data[key] = data[key]
  }
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
