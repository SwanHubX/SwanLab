<template>
  <div class="w-full">
    <EnvItems :data="item" v-for="item in environments" :key="item" />
    <EnvgpuItem />
  </div>
</template>

<script setup>
/**
 * @description: 系统硬件信息
 * @file: EnvHardware.vue
 * @since: 2024-01-24 21:19:24
 * 第一步优化已经做好了，怎么展示多个gpu信息，并且分别展示memory和tpye还需要再看
 **/

import { computed } from 'vue'
import { useExperimentStroe } from '@swanlab-vue/store'
import EnvItems from '../components/EnvItems.vue'
import EnvgpuItem from '../components/EnvgpuItem.vue'
const experimentStore = useExperimentStroe()
const experiment = experimentStore.experiment
const system = experiment.system
const environments = computed(() => {
  return [cpu.value, gpu.value]
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
      value: system.memory ? system.memory.toFixed(2) + 'GB' : ''
    }
  ]
})

const gpu = computed(() => {
  return [
    {
      key: 'gpu_cores',
      value: system.gpu?.cores
    },
    {
      key: 'gpu_type',
      value: system.gpu?.type[0]
    }
  ]
})
</script>

<style lang="scss" scoped></style>
