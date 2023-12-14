<template>
  <ErrorLayout>
    <component :is="InitError"></component>
  </ErrorLayout>
</template>

<script setup>
/**
 * @description: 错误页
 * @file: ErrorView.vue
 * @since: 2023-12-14 14:49:58
 **/
import ErrorLayout from '@swanlab-vue/layouts/ErrorLayout.vue'
import InitError from './pages/InitError.vue'

defineProps({
  code: {
    type: Number,
    default: 404
  }
})
</script>

<style lang="scss" scoped></style>
