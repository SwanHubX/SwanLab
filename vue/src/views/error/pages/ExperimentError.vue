<template>
  <ErrorWrap>
    <template #title>
      {{ $t(`error.experiment-error.title`) }}
    </template>
    <ul>
      <li v-for="index in [0, 1]" :key="index">
        {{ $t(`error.experiment-error.reasons[${index}]`) }}
      </li>
    </ul>
  </ErrorWrap>
</template>

<script setup>
/**
 * @description: 实验加载失败
 * @file: ExperimentError.vue
 * @since: 2023-12-15 15:18:59
 **/
import ErrorWrap from '../components/ErrorWrap.vue'
</script>

<style lang="scss" scoped></style>
