<template>
  <!-- 图表顶层，嵌入全局功能 -->
  <div></div>
  <!-- 每一个namespace对应一个图表容器 -->
  <ChartsContainer
    v-for="group in groups"
    :ref="(el) => setChartsRefList(el, index)"
    :key="group.name"
    :label="getNamespaces(group.name)"
    :charts="getCharts(group)"
  />
</template>

<script setup>
/**
 * @description: 图表面板，组织图表
 * @file: ChartsDashboard.vue
 * @since: 2024-03-04 19:30:29
 **/
import ChartsContainer from './ChartsContainer.vue'
import { t } from '@swanlab-vue/i18n'
import { ref } from 'vue'
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
</script>

<style lang="scss" scoped></style>
