<template>
  <div class="h-screen w-full flex justify-between">
    <!-- 左侧导航 -->
    <div class="w-80 bg-higher border-r flex-shrink-0 h-full flex flex-col">
      <slot name="left"></slot>
    </div>
    <!-- 右侧主体区域 -->
    <div class="w-full h-full overflow-y-auto">
      <slot></slot>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 默认布局，分为左测导航栏与右侧主体
 * @file: MainLayout.vue
 * @since: 2023-12-04 16:01:16
 **/
</script>

<style lang="scss" scoped></style>
