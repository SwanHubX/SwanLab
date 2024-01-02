<template>
  <div
    class="relative flex items-center justify-center rounded-full overflow-hidden animate-spin"
    :class="size ? `w-${size} h-${size}` : 'w-8 h-8'"
  >
    <div class="absolute w-full h-full rounded-full border-4 border-t-slate-500 border-slate-200"></div>
    <div class="w-3 h-3 z-10 rounded-full bg-slate-400"></div>
  </div>
</template>

<script setup>
/**
 * @description: 加载动画
 * @file: SLLoading.vue
 * @since: 2023-12-12 17:35:53
 **/

defineProps({
  size: {
    type: String,
    default: '8'
  }
})
</script>

<style lang="scss" scoped></style>
