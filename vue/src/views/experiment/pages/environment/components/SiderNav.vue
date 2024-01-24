<template>
  <div class="w-48">
    <router-link v-for="item in navs" :key="item.key" :to="item.to" active-class="active-router">
      {{ $t(`experiment.env.navs.${item.key}`) }}
    </router-link>
  </div>
</template>

<script setup>
/**
 * @description: 环境页 - 侧边导航栏
 * @file: SiderNav.vue
 * @since: 2024-01-09 15:32:18
 **/

const navs = [
  {
    key: 'index',
    to: `index`
  },
  {
    key: 'hardware',
    to: `hardware`
  },
  {
    key: 'requirements',
    to: `requirements`
  }
]
</script>

<style lang="scss" scoped>
a {
  @apply flex items-center px-4 h-11 text-sm text-default hover:bg-positive-dimmest rounded-lg mb-2;
}

.active-router {
  @apply bg-positive-dimmest text-positive-higher;
}
</style>
