<template>
  <div class="w-full">
    <EnvItems :data="item" v-for="item in environments" :key="item" />
  </div>
</template>

<script setup>
/**
 * @description: 实验环境 - 总览页
 * @file: IndexPage.vue
 * @since: 2024-01-09 15:47:20
 *
 * 这里会放一些实验配置信息
 **/

import { useExperimentStroe } from '@swanlab-vue/store'
import { computed } from 'vue'
import { formatTime } from '@swanlab-vue/utils/time'
import EnvItems from './EnvItems.vue'

const experimentStore = useExperimentStroe()
const experiment = experimentStore.experiment
const system = experiment.system

// ---------------------------------- 配置信息 ----------------------------------

const environments = computed(() => {
  return [times.value, systems.value]
})

// 时间相关
const times = computed(() => {
  return [
    {
      key: 'start_time',
      value: formatTime(experiment.create_time)
    },
    {
      key: 'duration',
      value: experimentStore.duration
    }
  ]
})

// 运行环境相关
const systems = computed(() => {
  return [
    {
      key: 'python',
      value: system.python
    },
    {
      key: 'executable',
      value: system.executable
    }
  ]
})
</script>

<style lang="scss" scoped></style>
