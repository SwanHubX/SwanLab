<template>
  <SLMenu class="right-0 w-96" corner-class="left-2.5">
    <button class="h-8 w-8 p-1 rounded text-dimmer hover:text-positive-default hover:bg-positive-dimmer">
      <SLIcon class="w-full h-full" icon="smooth" />
    </button>
    <template #pop> </template>
  </SLMenu>
</template>

<script setup>
/**
 * @description: 平滑按钮，包含菜单
 * @file: SmoothButton.vue
 * @since: 2024-03-04 21:13:50
 **/
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import SLMenu from '@swanlab-vue/components/menu/SLMenu.vue'
</script>

<style lang="scss" scoped></style>
