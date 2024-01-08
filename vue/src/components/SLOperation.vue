<template>
  <div class="absolute">
    <button class="outline-none" @click="show">
      <slot></slot>
    </button>
    <div
      class="absolute opacity-0 transition-all border p-2 rounded-lg"
      :class="`${isShow ? 'opacity-100' : ''} ${position}`"
    >
      <slot name="main"></slot>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 操作集合
 * @file: SLOperation.vue
 * @since: 2024-01-08 13:32:03
 *
 * 当有多个按钮、内容需要集成，并且具有收起、展开样式，可以使用该组件
 **/

import { ref } from 'vue'

defineProps({
  // 展开方向
  direction: {
    type: String,
    default: 'down',
    validator: (value) => {
      return ['up', 'down', 'right', 'left'].includes(value)
    }
  },
  // 相对位置
  position: {
    type: String,
    default: 'right-0 top-8'
  }
})

const isShow = ref(false)

const show = () => {
  isShow.value = !isShow.value
}
</script>

<style lang="scss" scoped></style>
