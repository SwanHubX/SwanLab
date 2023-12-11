<template>
  <ExtendBlock class="py-6" icon="experiment" :title="$t('experiment.index.config.title')">
    <div class="pl-6 w-full flex justify-between">
      <div class="w-1/2">
        <p class="font-semibold py-4">{{ $t('experiment.index.config.detail') }}</p>
        <div class="w-full">1</div>
      </div>
      <div class="w-1/2">
        <p class="font-semibold py-4">{{ $t('experiment.index.config.summarize') }}</p>
        <div class="w-full">2</div>
      </div>
    </div>
  </ExtendBlock>
</template>

<script setup>
/**
 * @description: 实验的配置信息
 * @file: ExperimentConfig.vue
 * @since: 2023-12-11 17:07:31
 **/
import ExtendBlock from '@swanlab-vue/views/experiment/components/ExtendBlock.vue'
</script>

<style lang="scss" scoped></style>
