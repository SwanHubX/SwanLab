<template>
  <ChartsPage :groups="groups" :charts="charts" v-if="groups.length" />
</template>

<script setup>
/**
 * @description: 项目对比图表，本组件完成项目对比图表的数据的请求和展示，大致流程是：
 * 1. 通过 http.get('/project/charts') 请求项目对比图表的数据，渲染到页面上
 * 2. 根据每个图表的数据源
 * @file: ChartsView.vue
 * @since: 2024-01-27 13:05:27
 **/
import http from '@swanlab-vue/api/http'
import { useProjectStore } from '@swanlab-vue/store'
import { ref, provide } from 'vue'
import ChartsPage from './components/ChartsPage.vue'
const projectStore = useProjectStore()
http.get('/project/charts').then(({ data }) => {
  // 将namespaces转换为groups
  charts.value = data.charts
  groups.value = data.namespaces.map((namespace) => {
    namespace.charts = namespace.charts.map((chart_id) => {
      // 在data.charts寻找id为当前chart_id的chart
      return data.charts.find((chart) => {
        return chart.id === chart_id
      })
    })
    return namespace
  })
})

// ---------------------------------- 数据驱动 ----------------------------------
// 项目对比图表数据，[{name, charts: [charts]}]
const groups = ref([])
const charts = ref([])

// ---------------------------------- 订阅管理器 ----------------------------------

class SubscriptionManager {
  constructor() {
    // 订阅管理器，map(experiment-id, Map(cid, callback))
    this.subscriptions = new Map()
  }
}

// ---------------------------------- 订阅函数 ----------------------------------

// ---------------------------------- 色盘注入 ----------------------------------
const colors = projectStore.colors
provide('colors', colors)
</script>

<style lang="scss" scoped></style>
