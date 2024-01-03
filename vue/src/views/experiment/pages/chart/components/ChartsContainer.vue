<template>
  <!-- 容器 -->
  <div class="charts-container">
    <!-- 控制栏 -->
    <div class="flex items-center">
      <!-- 关闭、展开按钮 -->
      <button class="flex items-center gap-1 relative" @click="handleExpand">
        <SLIcon icon="down" class="w-6 h-6 absolute -left-6" :class="{ '-rotate-90': !isExpand }" />
        {{ label }}
        <span class="px-3 bg-positive-dimmest text-sm rounded-full ml-2 text-positive-higher">{{ count || 0 }}</span>
      </button>
      <!-- 中间其他操作区 -->
      <div class="grow"></div>
      <!-- 添加图表的按钮 -->
    </div>
    <!-- 图表插槽 -->
    <div class="grid grid-cols-1 lg:grid-cols-2 xl:grid-cols-3 gap-6 mt-4.5" v-show="isExpand">
      <slot :isExpand="isExpand"></slot>
    </div>
  </div>
</template>

<script setup>
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import { ref } from 'vue'
/**
 * @description: 图表容器，用于包裹 PlotlyChart 组件，插槽允许isExpand参数，用于标识图表的展开和关闭
 * @file: ChartsContainer.vue
 * @since: 2023-12-09 15:35:28
 **/
defineProps({
  label: {
    type: String,
    required: true
  },
  count: {
    type: Number
  }
})

// ---------------------------------- 控制展开和关闭的状态 ----------------------------------

const isExpand = ref(true)

const handleExpand = () => {
  isExpand.value = !isExpand.value
}
</script>

<style lang="scss" scoped>
.charts-container {
  @apply border-b py-6 px-8;
}
</style>
