<template>
  <div class="group relative">
    <div class="absolute top-[-30px] left-[40%] text-center hidden group-hover:block">
      <p class="border px-3 py-1 text-sx scale-75 bg-higher rounded-lg">{{ `${line.key}: ${line.value}` }}</p>
    </div>
    <slot></slot>
  </div>
</template>

<script setup>
/**
 * @description: 表格行
 * @file: TableLine.vue
 * @since: 2024-01-24 16:40:01
 **/

defineProps({
  line: {
    type: Object,
    default: () => {
      return {}
    }
  }
})
</script>

<style lang="scss" scoped></style>
