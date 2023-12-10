<template>
  <!-- 标题 -->
  <div class="flex items-center gap-2 py-5 px-8 border-b">
    <SLIcon icon="experiment" class="w-5 h-5" />
    <h1 class="text-lg font-semibold">{{ projectStore.name + '/' + experiment.name }}</h1>
    <StatusLabel
      :id="experimentId"
      :status="experimentStatus"
      :name="experiment.name"
      v-if="experimentStatus !== undefined"
    />
  </div>
  <!-- 图表容器 -->
  <ChartsContainer label="default" :key="experimentId">
    <PlotlyChart v-for="(tag, index) in tags" :key="index" :sources="[tag]" />
  </ChartsContainer>
</template>

<script setup>
/**
 * @description: 实验-图表页
 * @file: ChartPage.vue
 * @since: 2023-12-09 20:39:41
 **/
import { inject, ref } from 'vue'
import { onBeforeRouteUpdate } from 'vue-router'
import { useProjectStore } from '@swanlab-vue/store'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import StatusLabel from '@swanlab-vue/components/StatusLabel.vue'
import ChartsContainer from './components/ChartsContainer.vue'
import PlotlyChart from './components/PlotlyChart.vue'
import http from '@swanlab-vue/api/http'

const projectStore = useProjectStore()

const experiment = inject('experiment')
const experimentId = inject('experimentId')

// ---------------------------------- 在此处完成路由监听和页面重新渲染 ----------------------------------
onBeforeRouteUpdate((to, from, next) => {
  console.log('刷新实验视图')
  next()
  // 先关闭轮询，然后重新开启一个
  clearInterval(timer)
  if (to.name === from.name) getExperiment(to.params.experimentId)
})

// ---------------------------------- 初始化：获取图表配置 ----------------------------------
const tags = ref()
const experimentStatus = ref(undefined)
// 用于设置轮询
let timer = undefined

/**
 * 根据id获取实验，如果id不传，默认使用experimentId
 * 如果实验正在进行，开启轮询
 * @param { string } id 实验id
 */
const getExperiment = (id = experimentId.value) => {
  http.get('/experiment/' + id).then(({ data }) => {
    // 判断data.tags和tags是否相同，如果不同，重新渲染
    if (data.tags.join('') !== tags.value?.join('')) {
      console.log('不同，重新渲染', data.tags, tags.value)
      tags.value = data.tags
    }
    // 如果status为0，且timer为undefined，开启轮询
    // 如果不为0且timer不为undefined，关闭轮询
    experimentStatus.value = data.status
    if (data.status === 0 && !timer) {
      timer = setInterval(() => {
        getExperiment(id)
      }, 5000)
    } else if (data.status !== 0 && timer) {
      clearInterval(timer)
    }
  })
}
// 立即执行
getExperiment()
</script>

<style lang="scss" scoped></style>
