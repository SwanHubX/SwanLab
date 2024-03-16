<template>
  <ChartsDashboard
    :groups="groups"
    :update-chart-status="updateChartStatus"
    :update-namespace-status="updateNamespaceStatus"
    :media="media"
    :default-color="defaultColor"
    :get-color="getColor"
    :subscribe="on"
  />
</template>

<script setup>
/**
 * @description: 项目对比图表，本组件完成项目对比图表的数据的请求和展示
 * 这个组件是为了实现自动刷新的无奈之举（在websocket没有上之前），整个组件刷新的时候，会重新请求数据
 * @file: ChartsPage.vue
 * @since: 2024-02-06 17:17:55
 **/
import { useProjectStore } from '@swanlab-vue/store'
import { onUnmounted } from 'vue'
import http from '@swanlab-vue/api/http'
import ChartsDashboard from '@swanlab-vue/charts/ChartsDashboard.vue'
import { updateChartStatus, updateNamespaceStatus, media } from '@swanlab-vue/api/chart'
const props = defineProps({
  // 图表组
  groups: {
    type: Array,
    required: true
  },
  // 图表列表
  charts: {
    type: Array,
    required: true
  },
  // 默认颜色
  defaultColor: {
    type: String,
    required: true
  },
  // 获取颜色
  getColor: {
    type: Function,
    required: true
  }
})

const projectStore = useProjectStore()

// ---------------------------------- 轮询器 ----------------------------------
const intervalMap = new Map()
const createInterval = (exp_name, cid, callback) => {
  const interval = setInterval(() => {
    getTagDataByExpName(exp_name, cid, callback)
  }, 5000)
  intervalMap.set(exp_name + '-' + cid, interval)
}

onUnmounted(() => {
  intervalMap.forEach((interval) => {
    clearInterval(interval)
  })
})

// ---------------------------------- 数据驱动 ----------------------------------

const on = (sources, cid, callback) => {
  // 获取数据
  sources.map((exp_name) => {
    // 如果exp_name对应的实验的show为0，则跳过
    if (projectStore.experiments.find((exp) => exp.name === exp_name).show === 0) return callback(exp_name, null, null)
    getTagDataByExpName(exp_name, cid, callback)
    // 判断当前实验状态
    if (projectStore.experiments.find((exp) => exp.name === exp_name).status === 0) {
      // 开启轮询
      createInterval(exp_name, cid, callback)
    }
  })
}

// ---------------------------------- 请求数据 ----------------------------------

/**
 * @description: 根据实验名称获取标签数据
 * @param {string} exp_name 实验名称
 * @param {number} cid 图表id
 */
const getTagDataByExpName = (exp_name, cid, callback) => {
  const exp_id = projectStore.experiments.find((exp) => exp.name === exp_name).id
  const tag_name = props.charts.find((chart) => chart.id === cid).name
  http
    .get(`/experiment/${exp_id}/tag/${tag_name}`)
    .then((res) => {
      callback(exp_name, res.data, null)
    })
    .catch(() => {
      callback(exp_name, null, null)
    })
}
</script>

<style lang="scss" scoped></style>
