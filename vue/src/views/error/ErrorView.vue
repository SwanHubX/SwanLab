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
import ExperimentError from './pages/ExperimentError.vue'
import { computed } from 'vue'
import { provide } from 'vue'

const props = defineProps({
  code: {
    type: Number,
    default: 404
  }
})

provide('error_code', props.code)

const Error = computed(() => {
  console.log(props.code)
  switch (props.code) {
    case 3500:
      return InitError
    case 3404:
      return ExperimentError
    default:
      return InitError
  }
})
</script>

<style lang="scss" scoped></style>
