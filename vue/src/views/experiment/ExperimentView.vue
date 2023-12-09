<template>
  <ExperimentLayout>
    <template #tabs>
      <TabsHeader />
    </template>
    <router-view />
  </ExperimentLayout>
</template>

<script setup>
/**
 * @description: 实验视图，展示实验列表
 * @file: ExperimentView.vue
 * @since: 2023-12-04 19:07:53
 **/
import TabsHeader from './components/TabsHeader.vue'
import ExperimentLayout from '@swanlab-vue/layouts/ExperimentLayout.vue'
</script>

<style lang="scss" scoped></style>
