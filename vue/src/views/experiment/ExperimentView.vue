<template>
  <!-- 标题 -->
  <div class="flex items-center gap-2 py-5 px-8 border-b">
    <SLIcon icon="experiment" class="w-5 h-5" />
    <h1 class="text-lg font-semibold">{{ projectStore.name + '/' + experiment.name }}</h1>
    <StatusLabel :id="experimentId" :status="experiment.status" :name="experiment.name" />
  </div>
  <!-- 图表容器 -->
  <ChartsContainer label="default" :key="experimentId">
    <PlotlyChart v-for="(tag, index) in tags" :key="index" :sources="[tag]" />
  </ChartsContainer>
</template>

<script setup>
/**
 * @description: 实验视图，展示实验列表
 * @file: ExperimentView.vue
 * @since: 2023-12-04 19:07:53
 **/
import { computed, provide, ref } from 'vue'
import { onBeforeRouteUpdate, useRoute } from 'vue-router'
// import http from '@swanlab-vue/api/http'
import { useProjectStore } from '@swanlab-vue/store'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import StatusLabel from '@swanlab-vue/components/StatusLabel.vue'
import ChartsContainer from './components/ChartsContainer.vue'
import PlotlyChart from './components/PlotlyChart.vue'
import http from '@swanlab-vue/api/http'
const route = useRoute()

const projectStore = useProjectStore()

const experimentId = computed(() => Number(route.params.experimentId))

// ---------------------------------- 在此处完成路由监听和页面重新渲染 ----------------------------------
onBeforeRouteUpdate((to, from, next) => {
  console.log('刷新实验视图')
  next()
  // 先关闭轮询，然后重新开启一个
  clearInterval(timer)
  console.log(timer)
  getExperiment(to.params.experimentId)
})

// ---------------------------------- 获取当前实验的配置 ----------------------------------
const experiment = computed(() => {
  return projectStore.experiments?.find((item) => item.experiment_id === experimentId.value)
})

const experimentColor = computed(() => {
  return experiment.value.color
})

// ---------------------------------- 初始化：获取图表配置 ----------------------------------
const tags = ref()
const experimentStatus = ref()
// 用于设置轮询
let timer = undefined

const getExperiment = (id = experimentId.value) => {
  http.get('/experiment/' + id).then(({ data }) => {
    // 判断data.tags和tags是否相同，如果不同，重新渲染
    if (data.tags.join('') !== tags.value?.join('')) {
      console.log('不同，重新渲染', data.tags, tags.value)
      tags.value = data.tags
    }
    // 如果status为0，且timer为undefined，开启轮询
    // 如果不为0且timer不为undefined，关闭轮询
    const status = data.status
    experimentStatus.value = ref()
    if (status === 0 && !timer) {
      timer = setInterval(() => {
        getExperiment(id)
      }, 5000)
    } else if (status !== 0 && timer) {
      clearInterval(timer)
    }
  })
}
// 立即执行
getExperiment()

// ---------------------------------- 依赖注入 ----------------------------------

provide('experimentId', experimentId)
provide('experimentColor', experimentColor)
</script>

<style lang="scss" scoped></style>
