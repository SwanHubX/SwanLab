<template>
  <section class="chart-container" @mouseenter="handleMouseEnter" @mouseleave="handleMouseLeave">
    <p class="text-center">{{ title }}</p>
    <div class="top-0 absolute" v-if="hover">dsadasdas</div>
    <slot />
  </section>
</template>

<script setup>
/**
 * @description: 图表容器，用于包裹图表组件，并且提供一些通用接口配置
 * @file: ChartContainer.vue
 * @since: 2023-12-12 21:01:26
 **/
defineProps({
  title: {
    type: String,
    required: true
  }
})

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
  flex: 1 0 400px;
}
</style>
