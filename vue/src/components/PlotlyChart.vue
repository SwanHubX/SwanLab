<template>
  <div class="w-full" ref="gd"></div>
</template>

<script setup>
/**
 * @description: 基于 Plotly.js 的图表组件，在此组件中封装
 * @file: PlotlyChart.vue
 * @since: 2023-12-09 00:13:00
 **/

import Plotly from 'plotly.js-dist-min'
import { onMounted } from 'vue'

const props = defineProps({
  data: {
    type: Object,
    required: true
  }
})

const gd = ref()

onMounted(() => {
  Plotly.newPlot(
    gd.value,
    /* JSON object */ {
      data: [{ y: props.data }]
    }
  )
})
</script>

<style lang="scss" scoped></style>
