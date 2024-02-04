<template>
  <div class="w-full border p-6 rounded-lg bg-default">
    <h1 class="w-full text-xl font-semibold pb-4 border-b mb-2">{{ $t(`experiment.env.title.${route.name}`) }}</h1>
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

import { useRoute } from 'vue-router'
import { useExperimentStore, useProjectStore } from '@swanlab-vue/store'
import { computed } from 'vue'
import { formatTime } from '@swanlab-vue/utils/time'
import EnvItems from '../components/EnvItems.vue'

const route = useRoute()
const experimentStore = useExperimentStore()
const projectStore = useProjectStore()
const experiment = experimentStore.experiment
const system = experiment.system

// ---------------------------------- 配置信息 ----------------------------------

// 环境配置汇总
const environments = computed(() => {
  return [overview.value, systems.value, gits.value, swanlab.value]
})

// 实验相关（实验描述、时间、系统）
const overview = computed(() => {
  return [
    // {
    //   key: 'name',
    //   value: experiment.name
    // },
    // {
    //   key: 'description',
    //   value: experiment.description
    // },
    {
      key: 'start_time',
      value: formatTime(experiment.create_time)
    },
    {
      key: 'duration',
      value: experimentStore.duration
    },
    {
      key: 'hostname',
      value: system.hostname
    },
    {
      key: 'OS',
      value: system.os
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
      highLight: true,
      copy: true
    },
    {
      key: 'workspace',
      value: system.cwd,
      highLight: true,
      copy: true
    },
    {
      key: 'command',
      value: system.command,
      highLight: true,
      copy: true
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
      value: system.git_info ? system.git_info[1] : '',
      copy: true
    }
  ]
})

// 硬件相关--v0.2.0 转移到System Hardware页面

// swanlab 相关
const swanlab = computed(() => {
  return [
    {
      key: 'logdir',
      value: projectStore.logdir + (projectStore.logdir.includes('/') ? '/' : '\\') + experiment.run_id,
      highLight: true,
      copy: true
    },
    {
      key: 'swanlab',
      value: `v${experiment.version}`
    }
  ]
})
</script>

<style lang="scss" scoped></style>
