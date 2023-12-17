<template>
  <ExperimentLayout v-if="ready" :key="experimentStore.id">
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
import { useExperimentStroe } from '@swanlab-vue/store'
import { onBeforeRouteUpdate, useRoute } from 'vue-router'
import http from '@swanlab-vue/api/http'
import { computed } from 'vue'
import { inject } from 'vue'
const route = useRoute()
const experimentStore = useExperimentStroe()
// ---------------------------------- 请求实验信息 ----------------------------------
const ready = computed(() => {
  return experimentStore.id !== undefined
})
const show_error = inject('show_error')
const init = async (id = route.params.experimentId) => {
  http
    .get(`/experiment/${id}`)
    .then(({ data }) => {
      experimentStore.experiment = data
    })
    .catch((response) => {
      // console.error(response)
      show_error(response.data.code)
    })
}

init()

onBeforeRouteUpdate((to, from) => {
  console.log('leave')
  if (to.params.experimentId !== from.params.experimentId) {
    init(to.params.experimentId)
  }
})
</script>

<style lang="scss" scoped></style>
