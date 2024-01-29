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
import { ref, inject } from 'vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import * as UTILS from './utils'

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

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------

const error = ref(props.chart.error)

// ---------------------------------- 图表颜色配置 ----------------------------------
// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const colors = inject('colors')
if (!colors) throw new Error('colors is not defined, please provide colors in parent component')
// ---------------------------------- 组件渲染逻辑 ----------------------------------

// ---------------------------------- 数据格式化 ----------------------------------

// ---------------------------------- 渲染、重渲染功能 ----------------------------------

// 渲染
const render = (data) => {
  console.log(data)
}
// 重渲染
const change = (data) => {}

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
