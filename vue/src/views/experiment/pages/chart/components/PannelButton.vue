<template>
  <button class="pannel-button" v-tippy="{ content: tip }">
    <SLIcon class="h-full aspect-square" :icon="icon" />
  </button>
</template>

<script setup>
import SLIcon from '@swanlab-vue/components/SLIcon.vue'

/**
 * @description: chart的按钮样式
 * @file: PannelButton.vue
 * @since: 2023-12-12 21:43:02
 **/

defineProps({
  icon: {
    type: String,
    required: true
  },
  tip: {
    type: String,
    required: true
  }
})
</script>

<style lang="scss" scoped>
.pannel-button {
  @apply text-dimmer h-full p-0.5 rounded hover:bg-positive-dimmest hover:text-positive-highest;
}
</style>
