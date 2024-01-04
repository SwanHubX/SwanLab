<template>
  <ErrorWrap>
    <template #title> Something Went Wrong </template>
    <p>{{ message }}</p>
  </ErrorWrap>
</template>

<script setup>
/**
 * @description: 自定义错误，用于一些不是很严重的错误
 * @file: CustomError.vue
 * @since: 2024-01-04 17:52:02
 **/
import { inject } from 'vue'
import ErrorWrap from '../components/ErrorWrap.vue'

const message = inject('error_message')
</script>

<style lang="scss" scoped></style>
