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
    <TextModule class="text-table" :data="data" v-if="originalData" :tag="title" />
    <!-- 放大效果弹窗 -->
    <SLModal class="pb-4 overflow-hidden" max-w="-1" v-model="isZoom" @onExit="exitByEsc">
      <TextModule v-model="isDetailZoom" :data="data" :tag="title" modal />
    </SLModal>
  </template>
</template>

<script setup>
/**
 * @description: 文字图标
 * @file: TextChart.vue
 * @since: 2024-02-20 20:04:09
 **/
import { ref, computed } from 'vue'
import SLModal from '@swanlab-vue/components/SLModal.vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
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

const originalData = ref()
const source = ref(props.chart.source)
const data = computed(() => {
  if (!originalData.value) return []
  return originalData.value[source.value[0]]?.list
})

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------

const error = computed(() => {
  if (!props.chart.error) {
    return props.chart.error
  } else if (!originalData.value) {
    return 'No Data'
  }
  return false
})

// ---------------------------------- 渲染、重渲染功能 ----------------------------------

// 渲染
const render = (data) => {
  // 保存原始数据，其中主要数据结构：{tag: {list: [{data: 'xxx', more: {caption: 'xxx'}}]}}
  originalData.value = data
}

// 重渲染
const change = (data) => {
  // 将发生更新的 tag 数据保存到原始数据中
  for (let key in data) {
    // originalData.value[key] = data[key] => 这行代码触发不了 props 的响应式，而下面这行可以
    originalData.value[key] = { ...data[key] }
  }
}

// ---------------------------------- 放大功能 ----------------------------------

// 是否放大
const isZoom = ref(false)
// 是否展示数据详情弹窗
const isDetailZoom = ref(false)
// 放大数据
const zoom = () => {
  isZoom.value = true
}

// 通过按键退出放大弹窗
const exitByEsc = () => {
  // 如果放大弹窗中有数据详情弹窗，则不关闭放大弹窗，直接关闭数据详情弹窗（在 TextModel.vue 中关闭）
  if (isDetailZoom.value) return
  isZoom.value = false
}

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  render,
  change,
  zoom
})
</script>

<style lang="scss" scoped>
.text-table {
  @apply pt-2 -mx-3 w-[calc(100%+1.5rem)];
}
</style>
