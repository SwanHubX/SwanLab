<template>
  <ExperimentLayout :key="experimentStore.id" v-if="ready">
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
import { useExperimentStroe, useProjectStore } from '@swanlab-vue/store'
import { onBeforeRouteUpdate, useRoute } from 'vue-router'
import http from '@swanlab-vue/api/http'
import { computed } from 'vue'
import { inject } from 'vue'
const route = useRoute()
const experimentStore = useExperimentStroe()
const projectStore = useProjectStore()
// ---------------------------------- 请求实验信息 ----------------------------------
const ready = computed(() => {
  return experimentStore.id !== undefined
})
const showErrorView = inject('showErrorView')
const init = async (id = route.params.experimentId) => {
  experimentStore.experiment = undefined
  // 清空charts
  experimentStore.charts = undefined
  http
    .get(`/experiment/${id}`)
    .then(({ data }) => {
      experimentStore.experiment = data
      // 设置实验状态
      projectStore.setExperimentStatus(data.experiment_id, data.status)
      // 如果实验状态还在running，就轮询
      if (experimentStore.isRunning) polling()
    })
    .catch((response) => {
      // console.error(response)
      showErrorView(response.data.code)
    })
}

init()

// ---------------------------------- 控制页面切换，组件刷新 ----------------------------------

onBeforeRouteUpdate((to, from) => {
  // console.log('leave')
  if (to.params.experimentId !== from.params.experimentId) {
    init(to.params.experimentId)
    clearInterval(timer)
  }
})

// ---------------------------------- 轮询函数，获取实验状态 ----------------------------------
let timer = null
const polling = () => {
  timer = setInterval(() => {
    http.get(`/experiment/${experimentStore.id}/status`).then(({ data }) => {
      if (Number(experimentStore.id) !== Number(route.params.experimentId)) {
        clearInterval(timer)
        return console.log('stop, experiment id is not equal to route id')
      }
      // 设置实验状态
      experimentStore.setStatus(data.status)
      projectStore.setExperimentStatus(experimentStore.id, data.status)
      experimentStore.charts = data.charts

      if (!experimentStore.isRunning) {
        clearInterval(timer)
        return console.log('stop, experiment status is not 0')
      }
    })
  }, 1000)
}
</script>

<style lang="scss" scoped></style>
