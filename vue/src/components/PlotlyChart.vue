<template>
  <div class="chart-container">
    <h1 class="plotly-title absolute z-10">{{ title }}</h1>
    <div class="plotly-chart absolute z-0" ref="gd"></div>
  </div>
</template>

<script setup>
/**
 * @description: 基于 Plotly.js 的图表组件，在此组件中封装
 * @file: PlotlyChart.vue
 * @since: 2023-12-09 00:13:00
 **/

import Plotly from 'plotly.js-dist-min'
import { onMounted, ref } from 'vue'

const props = defineProps({
  data: {
    type: Object,
    required: true
  },
  color: {
    type: String,
    required: true
  },
  title: {
    type: String,
    required: true
  }
})

const gd = ref()

onMounted(() => {
  Plotly.newPlot(
    gd.value,
    [{ y: props.data, line: { color: props.color } }],
    { showlegend: false },
    { displayModeBar: false, responsive: true }
  )
})
</script>

<style lang="scss" scoped>
.plotly-chart {
  height: calc(100% + 100px);
  width: calc(100% + 100px);
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
}
.plotly-title {
  top: 10%;
  left: 50%;
  transform: translate(-50%, -50%);
}

.chart-container {
  @apply w-full h-72 border rounded relative overflow-hidden;
  flex: 1 0 400px;
}
</style>
