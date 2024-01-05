<template>
  <div class="w-full flex flex-col items-center pt-10">
    <h1 class="text-5xl font-semibold">404</h1>
    <p class="py-4 text-dimmer font-semibold text-lg">{{ $t('error.404.message') }}</p>
    <SLButton class="px-4 py-1.5 border rounded-lg" @click="$router.replace('/')" theme="positive">{{
      $t('error.404.return')
    }}</SLButton>
  </div>
</template>

<script setup>
/**
 * @description: 路由不存在，对应路由 404
 * @file: 404.vue
 * @since: 2023-12-15 15:47:55
 **/

import SLButton from '@swanlab-vue/components/SLButton.vue'
</script>

<style lang="scss" scoped></style>
