<template>
  <div class="w-80 h-80" ref="chartRef"></div>
</template>

<script setup>
/**
 * @description: 测试ECharts的使用
 * @file: EChartsTest.vue
 * @since: 2023-11-24 21:40:59
 **/
import { onMounted, ref, watch } from 'vue'
import * as echarts from 'echarts'
const chartRef = ref()

const props = defineProps({
  data: {
    type: Array,
    default: () => []
  }
})

// 需要在onMounted中初始化图表，否则会报错
onMounted(() => {
  const chart = echarts.init(chartRef.value)
  chart.setOption({
    title: {
      text: 'ECharts 入门示例'
    },
    xAxis: {
      type: 'category',
      data: ['1', '2', '3']
    },
    yAxis: {
      type: 'value'
    },
    series: [
      {
        data: props.data,
        type: 'line'
      }
    ]
  })

  // 监听props.data的变化，一旦变化，就更新图表
  watch(
    () => props.data,
    (newVal) => {
      console.log('新的数据：', newVal)
      // 生成x轴的数据，0到newVal.length-1
      const xData = Array.from({ length: newVal.length }, (v, k) => k)
      console.log('新的x轴数据：', xData)
      chart.setOption({
        xAxis: {
          type: 'category',
          data: xData
        },
        series: [
          {
            data: newVal,
            type: 'line'
          }
        ]
      })
    }
  )
})
</script>

<style lang="scss" scoped></style>
