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
import EnvItems from '../components/EnvItems.vue'

const experimentStore = useExperimentStroe()
const experiment = experimentStore.experiment
const system = experiment.system

// ---------------------------------- 配置信息 ----------------------------------

// 环境配置汇总
const environments = computed(() => {
  return [times.value, systems.value, gits.value, hardware.value, swanlab.value]
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

// 系统环境相关
const systems = computed(() => {
  return [
    {
      key: 'python',
      value: system.python
    },
    {
      key: 'executable',
      value: system.executable,
      highLight: true
    },
    {
      key: 'command',
      value: system.command,
      highLight: true
    }
  ]
})

// git 相关
const gits = computed(() => {
  return [
    {
      key: 'git_remote',
      value: system.git_remote,
      highLight: true,
      link: true
    },
    {
      key: 'git_branch',
      value: system.git_info ? system.git_info[0] : ''
    },
    {
      key: 'git_commit',
      value: system.git_info ? system.git_info[1] : ''
    }
  ]
})

// 硬件相关
const hardware = computed(() => {
  return [
    {
      key: 'hostname',
      value: system.hostname
    },
    {
      key: 'OS',
      value: system.os
    },
    {
      key: 'memory',
      value: system.memory ? system.memory.toFixed(2) + 'GB' : ''
    },
    {
      key: 'cpu',
      value: system.cpu
    },
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

// swanlab 相关
const swanlab = computed(() => {
  return [
    {
      key: 'swanlab',
      value: `v${experiment.version}`
    }
  ]
})
</script>

<style lang="scss" scoped></style>
