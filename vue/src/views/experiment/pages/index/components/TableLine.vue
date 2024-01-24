<template>
  <div class="group w-full relative hover:bg-highest">
    <p
      class="border absolute max-w-[80%] left-10 translate-y-[-90%] top-0 text-center group-hover:block hidden px-3 py-1 text-sx scale-75 bg-higher rounded-lg"
    >
      {{ `${line.key}: ${line.value}` }}
    </p>
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
