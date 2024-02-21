<template>
  <div class="w-full h-full">
    <div class="w-full border" v-for="item in data" :key="item">
      {{ item }}
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 文字表格，被TextChart.vue调用
 * @file: TextModule.vue
 * @since: 2024-02-20 20:06:45
 **/

const props = defineProps({
  data: {
    type: Object,
    default: () => {}
  },
  source: {
    type: Array,
    default: () => []
  }
})
</script>

<style lang="scss" scoped></style>
