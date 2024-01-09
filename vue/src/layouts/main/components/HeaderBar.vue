<template>
  <div class="w-full h-full flex justify-between items-center px-6">
    <div class="flex items-center">
      <SLIcon class="w-8 h-8" icon="logo"></SLIcon>
      <span class="font-semibold text-2xl leading-6 mr-0.5">Swanlab</span>
      <span class="whitespace-nowrap text-xs">{{ version }}</span>
    </div>
  </div>
</template>

<script setup>
/**
 * @description:
 * @file: PageHeader.vue
 * @since: 2024-01-09 11:13:20
 **/

import SLIcon from '@swanlab-vue/components/SLIcon.vue'

defineProps({
  version: {
    type: String,
    default: 'unknown'
  }
})
</script>

<style lang="scss" scoped></style>
