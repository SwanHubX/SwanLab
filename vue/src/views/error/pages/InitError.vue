<template>
  <ErrorWrap>
    <template #title>
      {{ $t(`error.init-error.title`) }}
    </template>
    <ul>
      <li v-for="index in [0, 1, 2]" :key="index">
        {{ $t(`error.init-error.reasons[${index}]`) }}
      </li>
    </ul>
  </ErrorWrap>
</template>

<script setup>
/**
 * @description:
 * @file: InitError.vue
 * @since: 2023-12-14 14:49:52
 **/
import ErrorWrap from '../components/ErrorWrap.vue'
</script>

<style lang="scss" scoped></style>
