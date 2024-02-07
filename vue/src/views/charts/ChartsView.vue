<template>
  <!-- 第一行内容，项目标题、实验标题、编辑按钮、删除按钮 -->
  <div class="project-title transition-marging duration-300 mt-5 pl-6" :class="{ 'ml-8': !isSideBarShow }">
    <div class="flex items-center gap-3">
      <!-- 项目标题/实验标题 -->
      <h1 class="text-2xl items-center gap-1 font-semibold max-w-md truncate">
        {{ projectStore.name }}
      </h1>
    </div>
  </div>
  <ChartsPage :groups="groups" :charts="charts" :key="chartsPageKey" v-if="groups.length" />
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
import { ref, provide, inject } from 'vue'
import ChartsPage from './components/ChartsPage.vue'
import { onUnmounted } from 'vue'
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
const isSideBarShow = inject('isSideBarShow')

// ---------------------------------- 数据驱动 ----------------------------------
// 项目对比图表数据，[{name, charts: [charts]}]
const groups = ref([])
const charts = ref([])
const chartsPageKey = ref(0)

// ---------------------------------- 向projectStore注册回调，当点击眼睛时执行此回调 ----------------------------------
const handleShowChange = () => {
  // 重新渲染页面
  chartsPageKey.value++
}

projectStore.registerChangeShowCallback(handleShowChange)

onUnmounted(() => {
  projectStore.destoryChangeShowCallback()
})

// ---------------------------------- 色盘注入 ----------------------------------
const createGetSeriesColor = () => {
  // 遍历所有实验，实验名称为key，实验颜色为value
  const colors = projectStore.colorMap
  return (exp_name) => {
    return colors[exp_name]
  }
}
const colors = [...projectStore.colors]
colors.getSeriesColor = createGetSeriesColor()
// console.log(colors)
provide('colors', colors)
</script>

<style lang="scss" scoped></style>
