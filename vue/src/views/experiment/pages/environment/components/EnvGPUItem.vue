<template>
  <div class="flex flex-col gap-y-4 pt-4">
    <!-- gpu-cores -->
    <div class="w-full flex">
      <span class="text-dimmer w-44">GPU count</span>
      <span>{{ gpu_cores }}</span>
    </div>
    <!-- gpu detail message -->
    <div class="w-full flex" v-for="gpu in gpuList" :key="gpu.key">
      <!-- gpu -->
      <span class="text-dimmer w-44">{{ gpu.key }}</span>
      <!-- type -->
      <span class="w-96 shrink">{{ gpu.type }}</span>
      <!-- memory -->
      <span class="w-50 shrink-0">{{ gpu.memory }}</span>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 一行分为三部分：gpu , gpu_type , gpu_memory，主要是展示作用
 * @file: EnvGPUItem.vue
 * @since: 2024-01-25 22:52:50
 **/

import { ref } from 'vue'

const gpu_cores = ref(4)

// 模拟数据
const gpuList = [
  { key: 'GPU0', type: 'NVIDIA GeForce RTX 3060 Laptop GPU', memory: '40.00GB' },
  { key: 'GPU1', type: 'NVIDIA GeForce GTX 1660 Ti', memory: '32.00GB' },
  { key: 'GPU2', type: 'AMD Radeon RX 6700 XT', memory: '16.00GB' },
  { key: 'GPU3', type: 'AMD Radeon RX 6700 XT', memory: '16.00GB' }
]
</script>

<style lang="scss" scoped></style>
