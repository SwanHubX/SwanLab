<template>
  <div class="w-full h-full flex px-6 py-6">
    <SiderNav class="shrink-0" />
    <div class="w-full pl-14">
      <h1 class="w-full text-xl font-semibold pb-4 border-b mb-7">123123</h1>
      <RouterView></RouterView>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 实验环境页
 * @file: EnvironmentPage.vue
 * @since: 2024-01-09 15:14:05
 **/

import SiderNav from './components/SiderNav.vue'
</script>

<style lang="scss" scoped>
a {
  @apply flex items-center px-4 h-11 text-sm text-default hover:bg-positive-dimmest rounded-lg mb-2;
}

.active-router {
  @apply bg-positive-dimmest text-positive-higher;
}
</style>
