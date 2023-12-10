<template>
  {{ experiment }}
</template>

<script setup>
/**
 * @description: 实验-概览页
 * @file: ExperimentIndex.vue
 * @since: 2023-12-09 20:37:45
 **/
import { inject } from 'vue'

const experiment = inject('experiment')
</script>

<style lang="scss" scoped></style>
