<template>
  <div class="flex flex-col gap-y-4 pt-4" v-if="system.gpu?.cores">
    <!-- gpu-cores -->
    <div class="w-full flex">
      <span class="text-dimmer w-44">GPU count</span>
      <span>{{ system.gpu?.cores }}</span>
    </div>
    <!-- gpu detail message -->
    <div class="w-full flex" v-for="(gpu, i) in system.gpu.memory" :key="i">
      <!-- gpu -->
      <span class="text-dimmer w-44">GPU{{ i }}</span>
      <!-- type -->
      <span class="w-96 break-words whitespace-wrap shrink">{{ system.gpu.type[i] }}</span>
      <!-- memory -->
      <span class="shrink-0 pl-8">{{ gpu }} GB</span>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 一行分为三部分：gpu , gpu_type , gpu_memory，主要是展示作用
 * @file: EnvGPUItem.vue
 * @since: 2024-01-25 22:52:50
 **/

import { useExperimentStore } from '@swanlab-vue/store'

const experimentStore = useExperimentStore()
const experiment = experimentStore.experiment
const system = experiment.system
</script>

<style lang="scss" scoped></style>
