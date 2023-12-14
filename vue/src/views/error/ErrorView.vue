<template>
  <ErrorLayout>
    <component :is="Error"></component>
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
import { computed } from 'vue'

const props = defineProps({
  code: {
    type: Number,
    default: 404
  }
})

const Error = computed(() => {
  switch (props.code) {
    case 404:
      return InitError
    default:
      return InitError
  }
})
</script>

<style lang="scss" scoped></style>
