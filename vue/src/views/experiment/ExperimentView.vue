<template>
  <div class="flex">
    <PlotlyChart :data="data" v-if="data" />
    <PlotlyChart :data="data" v-if="data" />
  </div>
</template>

<script setup>
/**
 * @description: 实验视图，展示实验列表
 * @file: ExperimentView.vue
 * @since: 2023-12-04 19:07:53
 **/
import PlotlyChart from '@swanlab-vue/components/PlotlyChart.vue'
import { ref, computed } from 'vue'
import { useRoute } from 'vue-router'
import http from '@swanlab-vue/api/http'

const route = useRoute()

const experimentId = computed(() => route.params.experimentId)
// ---------------------------------- 在此处完成路由监听和页面重新渲染 ----------------------------------

// ---------------------------------- 获取当前实验的颜色配置 ----------------------------------

// ---------------------------------- 先模拟在这请求数据 ----------------------------------

const data = ref()

http.get(`/experiment/${experimentId.value}/test`).then((res) => {
  data.value = res.data.list
})
</script>

<style lang="scss" scoped></style>
