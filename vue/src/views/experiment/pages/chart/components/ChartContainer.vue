<template>
  <section class="chart-container" @mouseenter="handleMouseEnter" @mouseleave="handleMouseLeave">
    <p class="text-center">{{ chart.source[0] }}</p>
    <div class="chart-pannel" v-if="hover">
      <slot name="pannel" />
    </div>
    <slot :hover="hover" />
  </section>
</template>

<script setup>
/**
 * @description: 图表容器，用于包裹并渲染图表组件
 * @file: ChartContainer.vue
 * @since: 2023-12-12 21:01:26
 **/
import { ref } from 'vue'
const props = defineProps({
  chart: {
    type: Object,
    required: true
  }
})

console.log(props.chart)

// ---------------------------------- 判断是否处于hover状态 ----------------------------------
const hover = ref(false)
const handleMouseEnter = () => (hover.value = true)
const handleMouseLeave = () => (hover.value = false)
</script>

<style lang="scss" scoped>
.chart-container {
  @apply w-full h-72 border rounded relative overflow-hidden;
  @apply px-3 py-4;
  @apply flex-col flex justify-between;
}

.chart-pannel {
  @apply absolute top-1 right-2   h-4;
  @apply flex justify-end gap-2;
}
</style>
