<template>
  <div class="w-full flex justify-between items-center">
    <SLSearch class="max-w-96" :placeholder="$t('experiment.index.header.table-bar.summary')"></SLSearch>
    <div class="flex gap-5 items-center">
      <SLCheck :label="$t('experiment.index.header.table-bar.summary')" />
      <SLButton theme="default" class="px-3 py-2 rounded-lg flex items-center gap-2">
        <SLIcon icon="home" class="w-4 h-4" />
        <span>{{ $t('experiment.index.header.table-bar.export') }}</span>
      </SLButton>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 表格上方功能栏：搜索、筛选、导出
 * @file: TableBar.vue
 * @since: 2024-01-16 17:52:07
 **/
import SLCheck from '@swanlab-vue/components/SLCheck.vue'
</script>

<style lang="scss" scoped></style>
