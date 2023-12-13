<template>
  <div class="w-full">
    <!-- 展开按钮和标题 -->
    <button class="flex items-center py-4" @click="() => (isExtend = !isExtend)">
      <SLIcon icon="down" class="w-5 h-5 mx-1 text-dimmest" :class="isExtend ? 'rotate-show' : 'rotate-hidden'" />
      <SLIcon :icon="icon" class="w-6 h-6 mr-2" />
      <span class="text-lg font-semibold">{{ title }}</span>
    </button>
    <!-- 内容 -->
    <div
      class="transition-all duration-200"
      :class="isExtend ? 'max-h-[500px] overflow-auto' : 'max-h-0 overflow-hidden'"
    >
      <slot></slot>
    </div>
  </div>
</template>

<script setup>
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import { ref } from 'vue'

/**
 * @description: 拓展块
 * @file: ExtendBlock.vue
 * @since: 2023-12-11 17:11:32
 **/
defineProps({
  title: {
    type: String,
    required: true
  },
  icon: {
    type: String,
    required: true
  }
})

const isExtend = ref(true)
</script>

<style lang="scss" scoped>
@keyframes hidden {
  0% {
    transform: rotate(0deg);
  }
  100% {
    transform: rotate(-90deg);
  }
}

@keyframes show {
  0% {
    transform: rotate(-90deg);
  }
  100% {
    transform: rotate(0deg);
  }
}

.rotate-hidden {
  animation: hidden 0.2s linear forwards;
}

.rotate-show {
  animation: show 0.2s linear forwards;
}
</style>
