<template>
  <div class="w-full h-full">
    <div class="w-full bg-higher h-16">
      <slot name="tabs"></slot>
    </div>
    <div class="experiment-content">
      <slot></slot>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 实验页布局（不含左侧侧边栏）
 * @file: ExperimentLayout.vue
 * @since: 2023-12-09 20:22:32
 **/
</script>

<style scoped lang="scss">
.experiment-content {
  @apply w-full overflow-y-auto;
  height: calc(100vh - 4rem);
}
</style>
