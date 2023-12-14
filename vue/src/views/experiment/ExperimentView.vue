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
import { computed, provide } from 'vue'
import { useProjectStore, useExperimentStroe } from '@swanlab-vue/store'
import { useRoute } from 'vue-router'
import http from '@swanlab-vue/api/http'
const route = useRoute()
const projectStore = useProjectStore()
const experimentStore = useExperimentStroe()

// ---------------------------------- 请求实验信息 ----------------------------------

;(async (id = route.params.experimentId) => {
  http.get(`/experiment/${id}`).then(({ data }) => {
    experimentStore.experiment = data
  })
})()

// ---------------------------------- 获取当前实验的配置 ----------------------------------
const experiment = computed(() => {
  return projectStore.experiments?.find((item) => item.experiment_id === experimentId.value)
})

const experimentColor = computed(() => {
  return experiment.value.color
})

const experimentId = computed(() => Number(route.params.experimentId))

provide('experiment', experiment)
provide('experimentId', experimentId)
provide('experimentColor', experimentColor)
</script>

<style lang="scss" scoped></style>
