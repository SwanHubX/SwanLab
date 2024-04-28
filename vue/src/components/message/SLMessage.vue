<template>
  <div class="flex justify-center">
    <div class="message-container" :class="[message.type, theme === 'dark' ? 'dark-theme' : '']">
      <!-- 成功图标 -->
      <SLIcon v-if="message.type === 'success'" icon="success" class="text-white-default" />
      <!-- 错误图标 -->
      <SLIcon v-if="message.type === 'error'" icon="error" class="text-white-default" />
      <!-- 警告图标 -->
      <SLIcon v-if="message.type === 'warning'" icon="info" class="text-white-default" />
      <p class="text-nowrap">{{ message.text }}</p>
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
  },
  theme: {
    type: String,
    default: 'default'
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

<style lang="scss" scoped>
.dark-theme {
  background-color: #27282e !important;
  color: white !important;
  border-color: #616568 !important;
}
.message-container {
  @apply px-3 h-10 border rounded inline-flex items-center gap-2 shadow-md bg-default text-sm;
  svg {
    @apply w-5 h-5;
  }
}
.success {
  svg {
    @apply text-positive-default;
  }
}

.error {
  svg {
    @apply text-negative-default;
  }
}

.warning {
  svg {
    @apply text-warning-dimmer;
  }
}
</style>
