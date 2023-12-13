<template>
  <Tippy trigger="mouseenter" placement="right">
    <div
      class="w-5 h-5 flex items-center justify-center text-sm rounded-full cursor-default pointer-events-none bg-slate-400 text-white font-semibold"
      :class="iconClass"
    >
      ?
    </div>
    <template #content>
      <div :class="contentClass">{{ content }}</div>
      <slot></slot>
    </template>
  </Tippy>
</template>

<script setup>
/**
 * @description: 放置一些提示信息
 * @file: SLHelp.vue
 * @since: 2023-12-12 20:12:07
 *
 * 只用content传入提示不好改样式，要想自定义样式最好用插槽
 **/
import { Tippy } from 'vue-tippy'

defineProps({
  iconClass: {
    type: String,
    default: ''
  },
  content: {
    type: String,
    default: ''
  },
  contentClass: {
    type: String,
    default: ''
  }
})
</script>

<style lang="scss" scoped></style>
