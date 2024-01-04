<template>
  <div
    class="relative whitespace-nowrap w-full h-full"
    :class="`${head ? ' bg-[#f6f8fa] hover:bg-slate-200 py-4' : 'hover:bg-blue-100'} ${highLight ? 'bg-blue-50' : ''}`"
    :style="{ width: width + 'px' }"
  >
    <slot></slot>
    <!-- 拖拽点 -->
    <div class="drag" v-if="head && drag"></div>
  </div>
</template>

<script setup>
import { computed } from 'vue'

/**
 * @description: 用来规范单元格的样式
 * @file: TableCell.vue
 * @since: 2023-12-29 15:11:34
 **/
const props = defineProps({
  width: {
    type: String,
    defatult: 'auto'
  },
  head: {
    type: Boolean,
    defatult: false
  },
  drag: {
    type: Boolean,
    defatult: false
  },
  highLight: {
    type: Boolean,
    defatult: false
  }
})
</script>

<style lang="scss" scoped>
.drag {
  @apply w-1 h-full absolute right-0 top-0 hover:bg-positive-higher cursor-col-resize;
}
</style>
