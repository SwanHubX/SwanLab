<template>
  <div class="w-full flex flex-col items-center pt-10">
    <h1 class="text-5xl font-semibold">404</h1>
    <p class="py-4 text-dimmer">{{ $t('error.404.message') }}</p>
    <button class="px-4 py-1.5 border rounded-lg bg-positive-dimmer text-white" @click="$router.replace('/')">
      {{ $t('error.404.return') }}
    </button>
  </div>
</template>

<script setup>
/**
 * @description: 路由不存在
 * @file: 404.vue
 * @since: 2023-12-15 15:47:55
 **/
</script>

<style lang="scss" scoped></style>
