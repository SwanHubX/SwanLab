<template>
  <RouterLink class="flex gap-2 items-center w-full" :to="`/experiment/${id}`">
    <circle class="w-3 h-3 rounded-full flex-shrink-0" :style="{ backgroundColor: $props.color }" />
    <p class="truncate" :title="name">{{ name }}</p>
  </RouterLink>
</template>

<script setup>
import { RouterLink } from 'vue-router'

/**
 * @description: 展示实验名称和实验对应的颜色，并且可以完成实验跳转
 * @file: ExperimentName.vue
 * @since: 2023-12-08 20:16:08
 **/

defineProps({
  name: {
    type: String,
    required: true
  },
  color: {
    type: String,
    required: true
  },
  id: {
    required: true
  }
})
</script>

<style lang="scss" scoped></style>
