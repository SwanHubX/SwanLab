<template>
  <!-- 标题 -->
  <div class="flex items-center gap-2 py-5 px-8 border-b">
    <div class="w-5 h-5 rounded-full" :style="{ backgroundColor: experimentStore.color }"></div>
    <h1 class="text-lg font-semibold">{{ experimentStore.name }}</h1>
    <SLStatusLabel :status="experimentStore.status" />
  </div>
  <!-- 图表容器 -->
  <template v-if="status === 'success'">
    <ChartsContainer
      v-for="group in groups"
      :key="group.namespace"
      :label="getGroupName(group.namespace)"
      :count="group.charts.length"
    >
      <ChartContainer :chart="charts[cnMap.get(chartID)]" v-for="chartID in group.charts" :key="chartID" />
    </ChartsContainer>
  </template>
  <EmptyExperiment v-once v-else-if="status !== 'initing'" />
</template>

<script setup>
/**
 * @description: 实验图表页，在本处将完成图表布局的初始化，注入订阅依赖等内容
 * @file: ChartPage.vue
 * @since: 2023-12-25 15:34:51
 **/
import SLStatusLabel from '@swanlab-vue/components/SLStatusLabel.vue'
import { useExperimentStroe } from '@swanlab-vue/store'
import http from '@swanlab-vue/api/http'
import { ref } from 'vue'
import ChartsContainer from './components/ChartsContainer.vue'
import ChartContainer from './components/ChartContainer.vue'
import EmptyExperiment from './components/EmptyExperiment.vue'
import { t } from '@swanlab-vue/i18n'
const experimentStore = useExperimentStroe()

// ---------------------------------- 主函数：获取图表配置信息，渲染图表布局 ----------------------------------
// 基于返回的namespcaes和charts，生成一个映射关系,称之为cnMap
// key是chart_id, value是chart_id在charts中的索引
let cnMap = new Map()
let charts, groups
const status = ref('initing')
;(async function () {
  const { data } = await http.get(`/experiment/${experimentStore.id}/chart`)
  charts = data.charts || []
  groups = data.namespaces || []
  // 添加一个临时array，用于减少遍历次数
  const chartsIdArray = charts.map((chart) => chart.chart_id)
  // 标志，判断groups是否找到对应的chart
  let found = false
  // 遍历groups
  groups.forEach((group) => {
    found = false
    group.charts.forEach((id) => {
      // 遍历charts，如果id和charts.charts_id相同，cnMap添加一个key
      for (let i = 0; i < chartsIdArray.length; i++) {
        const chartId = chartsIdArray[i]
        if (chartId === id) {
          found = true
          cnMap.set(chartId, i)
          break
        }
      }
      if (!found) {
        console.error('charts: ', charts)
        console.error('ngroups: ', groups)
        throw new Error('Charts and groups cannot correspond.')
      }
    })
  })
  // 遍历完了，此时cnMap拿到了,模版部分可以依据这些东西渲染了
  // 但是还需要完成订阅模型的初始化
  // console.log(cnMap)
})().then(() => {
  status.value = 'success'
})

// ---------------------------------- 发布、订阅模型 ----------------------------------

// ---------------------------------- 工具函数：辅助渲染 ----------------------------------

/**
 * 将组名进行一些翻译
 * @param { string } namespace 组名
 */
const getGroupName = (namespace) => {
  console.log(namespace)
  if (namespace === 'default') return t('experiment.chart.label.default')
  else return namespace
}
</script>

<style lang="scss" scoped></style>
