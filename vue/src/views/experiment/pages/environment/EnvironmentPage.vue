<template>
  <div
    class="w-full h-full flex p-6 bg-higher"
    :class="experimentStore.experiment.description ? 'min-h-[calc(100vh-195px)]' : 'min-h-[calc(100vh-163px)]'"
  >
    <SiderNav class="shrink-0 mr-6" />
    <div class="w-full rounded-lg mr-6">
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
import { useExperimentStore } from '@swanlab-vue/store'

const experimentStore = useExperimentStore()
</script>

<style lang="scss" scoped></style>
