<template>
  <ChartsContainer v-for="group in groups" :key="group.name" :label="group.name" :charts="group.charts" />
</template>

<script setup>
/**
 * @description: 项目对比图表，本组件完成项目对比图表的数据的请求和展示
 * 这个组件是为了实现自动刷新的无奈之举（在websocket没有上之前），整个组件刷新的时候，会重新请求数据
 * @file: ChartsPage.vue
 * @since: 2024-02-06 17:17:55
 **/
import ChartsContainer from '@swanlab-vue/charts/ChartsContainer.vue'
import { useProjectStore } from '@swanlab-vue/store'
import { provide } from 'vue'
import http from '@swanlab-vue/api/http'
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
  }
})
const projectStore = useProjectStore()

// ---------------------------------- 数据驱动 ----------------------------------

provide('$on', (sources, cid, callback) => {
  // 获取数据
  sources.map((exp_name) => {
    const exp_id = projectStore.experiments.find((exp) => exp.name === exp_name).id
    const tag_name = props.charts.find((chart) => chart.id === cid).name
    return new Promise((resolve, reject) => {
      http
        .get(`/experiment/${exp_id}/tag/${tag_name}`)
        .then((res) => {
          callback(exp_name, res.data, null)
          resolve()
        })
        .catch((err) => {
          // callback(exp_name, null, err)
          callback(exp_name, null, null)
          resolve()
        })
    })
  })
})

provide('$off', (sources, cid, callback) => {})

// ---------------------------------- 请求数据 ----------------------------------

/**
 * 通过实验id获取标签数据，
 * @param {  } expId
 */
const getTagDataByExpId = (expId) => {
  return http.get(`/experiment/${expId}/tags`)
}
</script>

<style lang="scss" scoped></style>
