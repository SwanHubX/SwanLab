<template>
  <div class="h-screen w-full">
    <!-- 顶部 -->
    <div class="w-full h-[52px] bg-higher">
      <HeaderBar :version="version" />
    </div>
    <div class="flex justify-between h-[calc(100vh-52px)]">
      <!-- 左侧导航 -->
      <div class="w-80 border-r flex-shrink-0 h-full flex flex-col">
        <HomeSiderBar />
      </div>
      <!-- 右侧主体区域 -->
      <div class="w-full h-full overflow-y-auto">
        <slot></slot>
      </div>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 默认布局，分为左测导航栏与右侧主体
 * @file: MainLayout.vue
 * @since: 2023-12-04 16:01:16
 **/

import HomeSiderBar from './HomeSiderBar.vue'
import HeaderBar from './HeaderBar.vue'

defineProps({
  version: {
    type: String,
    default: ''
  }
})
</script>

<style lang="scss" scoped></style>
