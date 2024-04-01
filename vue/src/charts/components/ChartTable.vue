<template>
  <div class="w-full h-full">
    {{ header }}
    <hr />
    {{ data }}
  </div>
</template>

<script setup>
/**
 * @description: 表格
 * @file: ChartTable.vue
 * @since: 2024-04-01 17:51:44
 **/
const props = defineProps({
  header: {
    type: Array,
    default: () => []
  },
  data: {
    type: Array,
    default: () => []
  }
})
</script>

<style lang="scss" scoped></style>
