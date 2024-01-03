<template>
  <div class="px-3 h-10 border border-gray bg-default rounded">
    <div class="flex items-center h-full">
      <p>{{ message.text }}</p>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 消息提示组件，用于显示单个消息提示，由SLMessages组件调用，不应该单独使用
 * @file: SLMessage.vue
 * @property { string } text 消息提示的信息
 * @property { string } type 消息提示的类型，可选值为error、success、processing
 * @since: 2024-01-03 15:34:01
 **/
import { onMounted, onUnmounted } from 'vue'

const props = defineProps({
  message: {
    type: Object,
    required: true
  }
})

let timer = null
onMounted(() => {
  timer = setTimeout(() => {
    props.message.close()
  }, props.message.delay)
})

onUnmounted(() => {
  clearTimeout(timer)
})
</script>

<style lang="scss" scoped></style>
