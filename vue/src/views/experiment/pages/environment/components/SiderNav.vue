<template>
  <div class="w-48">
    <router-link v-for="item in navs" :key="item.key" :to="item.to" active-class="active-router" exact>
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

import { useExperimentStroe } from '@swanlab-vue/store'

const id = useExperimentStroe().id

const navs = [
  {
    key: 'index',
    to: `/experiment/${id}/env/index`
  },
  {
    key: 'dependances',
    to: `/experiment/${id}/env/dependencies`
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
