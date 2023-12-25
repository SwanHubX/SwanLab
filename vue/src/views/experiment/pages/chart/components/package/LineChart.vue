<template>
  <div ref="g2Ref"></div>
  <!-- 放大效果 -->
  <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom">
    <p class="text-center mt-4 mb-10 text-2xl font-semibold">{{ title }}</p>
    <div ref="g2ZoomRef"></div>
    <p class="absolute right-12 bottom-16 text-xs text-dimmer scale-90">{{ tagData.xTitle }}</p>
  </SLModal>
</template>

<script setup>
/**
 * @description: 折线图表
 * @file: LineChart.vue
 * @since: 2023-12-25 20:17:19
 **/
import SLModal from '@swanlab-vue/components/SLModal.vue'
import { ref } from 'vue'

// 是否放大
const isZoom = ref(false)
</script>

<style lang="scss" scoped></style>
