<template>
  <div class="flex flex-col justify-center grow text-dimmer gap-2">
    <SLIcon class="mx-auto h-5 w-5" icon="error" />
    <p class="text-center text-xs">{{ $t('chart.unkonwn', { type: chart.type }) }}</p>
  </div>
</template>

<script setup>
/**
 * @description: 未知图表，用于展示未知类型的图表
 * @file: UnknownChart.vue
 * @since: 2023-12-29 15:55:18
 **/

defineProps({
  title: {
    type: String,
    required: true
  },
  chart: {
    type: Object,
    required: true
  }
})

const render = () => {}
const change = () => {}
const zoom = () => {}

defineExpose({
  render,
  change,
  zoom
})
</script>

<style lang="scss" scoped></style>
