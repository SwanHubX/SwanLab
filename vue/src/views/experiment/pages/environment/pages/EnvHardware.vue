<template>
  <div class="w-full border p-6 rounded-lg bg-default">
    <h1 class="w-full text-xl font-semibold pb-4 border-b mb-2">{{ $t(`experiment.env.title.${route.name}`) }}</h1>
    <EnvItems :data="item" v-for="item in environments" :key="item" />
    <EnvGPUItem />
    <div v-if="Object.keys(experimentStore.experiment.system).length === 0">
      <p class="text-center pt-5">No hardware information</p>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 系统硬件信息
 * @file: EnvHardware.vue
 * @since: 2024-01-24 21:19:24
 **/

import { computed } from 'vue'
import { useExperimentStore } from '@swanlab-vue/store'
import EnvItems from '../components/EnvItems.vue'
import EnvGPUItem from '../components/EnvGPUItem.vue'
import { useRoute } from 'vue-router'
const experimentStore = useExperimentStore()
const experiment = experimentStore.experiment
const system = experiment.system
const route = useRoute()
const environments = computed(() => {
  return [cpu.value]
})

// 系统硬件信息
const cpu = computed(() => {
  return [
    {
      key: 'cpu',
      value: system.cpu
    },
    {
      key: 'memory',
      value: system.memory ? system.memory + 'GB' : ''
    }
  ]
})
</script>

<style lang="scss" scoped></style>
