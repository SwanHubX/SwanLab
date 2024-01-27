<template>
  <div>
    <router-link v-for="item in navs" :key="item.key" :to="item.to" active-class="active-router">
      <span class="lg:w-48 lg:h-auto lg:px-4 lg:transition-none px-6 w-0 h-0 transition-all overflow-hidden">{{
        $t(`experiment.env.navs.${item.key}`)
      }}</span>
      <SLIcon :icon="item.icon" class="w-4 h-4 lg:w-0 lg:h-0 transition-all absolute left-4" />
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
    to: `index`,
    icon: 'home'
  },
  {
    key: 'hardware',
    to: `hardware`,
    icon: 'config'
  },
  {
    key: 'requirements',
    to: `requirements`,
    icon: 'help'
  }
]
</script>

<style lang="scss" scoped>
a {
  @apply relative flex items-center h-11 text-default hover:bg-highest rounded-lg mb-2;
}

.active-router {
  @apply bg-highest hover:bg-highest font-semibold;
}
</style>
