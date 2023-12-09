<template>
  <!-- 标题 -->
  <div class="flex items-center gap-2 py-5 px-8 border-b">
    <SLIcon icon="experiment" class="w-5 h-5" />
    <h1 class="text-lg font-semibold">{{ projectStore.name + '/' + experiment.name }}</h1>
    <StatusLabel :id="experimentId" :status="experiment.status" :name="experiment.name" />
  </div>
  <!-- 图表容器 -->
  <ChartsContainer label="default">
    <PlotlyChart :data="data" :color="experiment.color" :title="tag" v-if="!!data" />
    <PlotlyChart :data="data" :color="experiment.color" :title="tag" v-if="!!data" />
    <PlotlyChart :data="data" :color="experiment.color" :title="tag" v-if="!!data" />
    <PlotlyChart :data="data" :color="experiment.color" :title="tag" v-if="!!data" />
    <PlotlyChart :data="data" :color="experiment.color" :title="tag" v-if="!!data" />
  </ChartsContainer>
</template>

<script setup>
/**
 * @description: 实验视图，展示实验列表
 * @file: ExperimentView.vue
 * @since: 2023-12-04 19:07:53
 **/
import PlotlyChart from '@swanlab-vue/components/PlotlyChart.vue'
import { ref, computed } from 'vue'
import { onBeforeRouteUpdate, useRoute } from 'vue-router'
import http from '@swanlab-vue/api/http'
import { useProjectStore } from '@swanlab-vue/store'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import StatusLabel from '@swanlab-vue/components/StatusLabel.vue'
import ChartsContainer from './components/ChartsContainer.vue'
const route = useRoute()

const projectStore = useProjectStore()

const experimentId = computed(() => Number(route.params.experimentId))
// ---------------------------------- 在此处完成路由监听和页面重新渲染 ----------------------------------
onBeforeRouteUpdate((to, from, next) => {
  console.log('刷新实验视图')
  next()
})

// ---------------------------------- 获取当前实验的配置 ----------------------------------
const experiment = computed(() => {
  return projectStore.experiments?.find((item) => item.experiment_id === experimentId.value)
})

// ---------------------------------- 先模拟在这请求数据 ----------------------------------

const data = ref()
const tag = 'test'
http.get(`/experiment/${experimentId.value}/` + tag).then((res) => {
  data.value = res.data.list
})
</script>

<style lang="scss" scoped></style>
