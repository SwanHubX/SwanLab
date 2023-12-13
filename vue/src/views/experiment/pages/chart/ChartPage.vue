<template>
  <!-- 标题 -->
  <div class="flex items-center gap-2 py-5 px-8 border-b">
    <SLIcon icon="experiment" class="w-5 h-5" />
    <h1 class="text-lg font-semibold">{{ experiment.name }}</h1>
    <SLStatusLabel :status="experimentStatus" v-if="experimentStatus !== undefined" />
  </div>
  <!-- 图表容器 -->
  <ChartsContainer label="default" :key="experimentId" v-if="showContainer(tags)">
    <!-- TODO 后续多数据源的时候，这里需要改变sources的保存逻辑 -->
    <G2Chart v-for="(tag, index) in tags" :key="index" :sources="[tag]" />
  </ChartsContainer>
</template>

<script setup>
/**
 * @description: 实验-图表页
 * @file: ChartPage.vue
 * @since: 2023-12-09 20:39:41
 **/
import { inject, ref, provide } from 'vue'
import { onBeforeRouteUpdate, onBeforeRouteLeave } from 'vue-router'
import SLStatusLabel from '@swanlab-vue/components/SLStatusLabel.vue'
import ChartsContainer from './components/ChartsContainer.vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import G2Chart from './components/G2Chart.vue'
import http from '@swanlab-vue/api/http'

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

onBeforeRouteLeave((to, from, next) => {
  console.log('离开实验视图')
  next()
  // 离开页面时，关闭轮询
  clearInterval(timer)
})

// ---------------------------------- 初始化：获取图表配置 ----------------------------------
const tags = ref()
const experimentStatus = ref(undefined)
provide('experimentStatus', experimentStatus)
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
      // console.log('不同，重新渲染', data.tags, tags.value)
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

// ---------------------------------- 判断某个namespcace下是否存在图表 ----------------------------------
const showContainer = (sources) => {
  return !!sources?.length
}
</script>

<style lang="scss" scoped></style>
