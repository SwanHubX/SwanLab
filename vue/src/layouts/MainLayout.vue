<template>
  <div class="min-h-screen w-screen flex justify-between">
    <!-- 左侧导航 -->
    <div class="min-w-[320px] bg-dimmer border-r border-default">
      <slot name="left"></slot>
    </div>
    <!-- 右侧主体区域 -->
    <div class="w-full h-full">
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
