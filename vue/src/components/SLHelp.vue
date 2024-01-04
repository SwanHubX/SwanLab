<template>
  <Tippy trigger="mouseenter" placement="right" :class="document ? 'cursor-pointer' : ''" @click="learn_more">
    <div
      class="w-5 h-5 flex items-center justify-center text-sm rounded-full pointer-events-none bg-overlay text-white-default font-semibold"
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

const props = defineProps({
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
  },
  document: {
    type: String
  }
})

const learn_more = () => {
  props.document ? window.open(props.document, '_blank') : ''
}
</script>

<style lang="scss" scoped></style>
