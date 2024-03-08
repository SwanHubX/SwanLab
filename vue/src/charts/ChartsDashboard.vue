<template>
  <!-- 图表顶层，嵌入全局功能 -->
  <div class="charts-pannel">
    <SmoothButton @smooth="handleSmooth" :default-method-id="1" />
  </div>
  <!-- 每一个namespace对应一个图表容器 -->
  <ChartsContainer
    v-for="(group, index) in groups"
    :ref="(el) => setChartsRefList(el, index)"
    :key="group.name"
    :label="getNamespaces(group.name)"
    :charts="getCharts(group)"
    :opened="!!group.opened"
    @switch="(opened) => debouncedHandleSwitch(group.id, opened)"
  />
</template>

<script setup>
/**
 * @description: 图表面板，组织图表
 * @file: ChartsDashboard.vue
 * @since: 2024-03-04 19:30:29
 **/
import ChartsContainer from './components/ChartsContainer.vue'
import { t } from '@swanlab-vue/i18n'
import SmoothButton from './components/SmoothButton.vue'
import { debounce } from '@swanlab-vue/utils/common'
import http from '@swanlab-vue/api/http'
const props = defineProps({
  // 整个图表列表集合
  groups: {
    type: Array,
    required: true
  }
})

const chartsRefList = ref([])

const setChartsRefList = (el, index) => {
  chartsRefList.value[index] = el
  chartsRefList.value.length = props.groups.length
}

// ---------------------------------- 获取namespace名称 ----------------------------------
/**
 * 将组名进行一些翻译
 * @param { string } namespace 组名
 */
const getNamespaces = (namespace) => {
  // console.log(name)
  if (namespace === 'default') return t('common.chart.label.default')
  else return namespace
}

// ---------------------------------- 获取图表列表 ----------------------------------
/**
 * 根据group中的chart_id获取chart配置，生成一个数组
 */
const getCharts = (group) => {
  return group.charts
}

// ---------------------------------- 处理图表平滑事件 ----------------------------------

const handleSmooth = (method) => {
  // 遍历chartsRefList中的所有chartRefList，调用其smooth方法
  chartsRefList.value.forEach((chartsRef) => {
    chartsRef.chartRefList.forEach((chartRef) => {
      chartRef.smooth(method)
    })
  })
}

// ---------------------------------- 在此处处理命名空间打开和关闭 ----------------------------------
const handleSwitch = (id, opened) => {
  console.log('namespace click', id, opened)
  // 向后端更新展开状态
  http.patch('/namespace/' + id + '/opened', {
    opened
  })
}

const debouncedHandleSwitch = debounce(handleSwitch, 300)
</script>

<style lang="scss" scoped>
.charts-pannel {
  @apply sticky px-4 py-3 bg-higher flex justify-end top-0 border-b;
  z-index: 30;
}
</style>
