<template>
  <!-- 图表顶层，嵌入全局功能 -->
  <div class="charts-pannel">
    <!-- 搜索框 -->
    <div class="w-60 mr-2">
      <SLSearch :placeholder="$t('chart.search.placeholder')" reverse @search="search" @input="search" />
    </div>
    <SmoothButton @smooth="handleSmooth" :default-method-id="1" />
  </div>
  <!-- 每一个namespace对应一个图表容器 -->
  <ChartsContainer
    v-for="(group, index) in nowChartsGroupList"
    :ref="(el) => setChartsRefList(el, index)"
    :key="group.id"
    :label="getNamespaces(group)"
    :charts="getCharts(group)"
    :index="index"
    :opened="!!group.opened"
    :isPinned="group.id === -1"
    :isHidden="group.id === -2"
    @switch="(opened) => debouncedHandleSwitch(opened, group)"
    @pin="pin"
    @unpin="unpin"
    @hide="hide"
    @unhide="unhide"
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
import { debounces } from '@swanlab-vue/utils/common'
import SLSearch from '@swanlab-vue/components/SLSearch.vue'
import { provide } from 'vue'

const props = defineProps({
  // 整个图表列表集合
  groups: {
    type: Array,
    required: true
  },
  // 更新图表状态
  updateChartStatus: {
    type: Function,
    required: true
  },
  // 更新命名空间状态
  updateNamespaceStatus: {
    type: Function,
    required: true
  },
  // 媒体获取对象
  media: {
    type: Object,
    required: true
  },
  // 图表默认颜色
  defaultColor: {
    type: String,
    required: true
  },
  // 获取某一个数据的颜色
  getColor: {
    type: Function,
    required: true
  },
  // 订阅数据更新
  subscribe: {
    type: Function,
    required: true
  },
  // 取消订阅数据更新
  unsubscribe: {
    type: Function,
    default: () => {}
  },
  // 图表跳转实验前缀
  experimentPrefix: {
    type: String,
    default: '/experiment/'
  }
})
const chartsGroupList = ref(props.groups)
watch(
  () => props.groups,
  (newVal) => {
    // 当groups发生变化时，更新chartsGroupList
    chartsGroupList.value = newVal
  }
)
const searchChartsGroupList = ref(null)

const chartsRefList = ref([])

const nowChartsGroupList = computed(() => {
  return searchChartsGroupList.value ? searchChartsGroupList.value : chartsGroupList.value
})

const setChartsRefList = (el, index) => {
  chartsRefList.value[index] = el
  chartsRefList.value.length = chartsGroupList.value.length
}

// ---------------------------------- 获取namespace名称 ----------------------------------
/**
 * 将组名进行一些翻译
 * @param { Object } namespace 命名空间配置
 */
const getNamespaces = (namespace) => {
  // console.log(name)
  if (namespace.name === 'default') return t('chart.label.default')
  if (namespace.id === -1) return t('chart.label.pinned')
  if (namespace.id === -2) return t('chart.label.hidden')
  else return namespace.name
}

// ---------------------------------- 获取图表列表 ----------------------------------
/**
 * 根据group中的chart_id获取chart配置，生成一个数组
 */
const getCharts = (group) => {
  return group.charts
}

// ---------------------------------- 处理图表平滑事件 ----------------------------------
const smoothMethod = ref(null)
const handleSmooth = (method) => {
  // 遍历chartsRefList中的所有chartRefList，调用其smooth方法
  chartsRefList.value.forEach((chartsRef) => {
    chartsRef.chartRefList.forEach((chartRef) => {
      chartRef?.smooth(method)
    })
  })
  smoothMethod.value = method
}
// 暂时先以provide的方式传递smoothMethod，后续会存储在后端，通过接口获取
provide('smoothMethod', smoothMethod)

// ---------------------------------- 在此处处理命名空间打开和关闭 ----------------------------------
const changeNamespaceStatus = (opened, namespace) => {
  // console.log('namespace click', opened)
  // 向后端更新展开状态
  props.updateNamespaceStatus(opened, namespace)
}

const debouncedHandleSwitch = debounces(changeNamespaceStatus, 300)

const changeChartStatus = async (chart, status) => {
  const data = await props.updateChartStatus(chart, status)
  chartsGroupList.value = data.groups
}

// ---------------------------------- 处理pin事件 ----------------------------------
/**
 * 图表置顶
 * @param { Object } chart 图表配置
 */
const pin = (chart) => {
  changeChartStatus(chart, 1)
}

/**
 * 图表取消置顶
 * @param { Object } chart 图表配置
 */
const unpin = (chart) => {
  changeChartStatus(chart, 0)
}

// ---------------------------------- 处理hide事件 ----------------------------------
/**
 * 图表隐藏
 * @param { Object } chart 图表配置
 */
const hide = (chart) => {
  changeChartStatus(chart, -1)
}

/**
 * 图表取消隐藏
 * @param { Object } chart 图表配置
 */
const unhide = (chart) => {
  changeChartStatus(chart, 0)
}

// ---------------------------------- 处理搜索图表事件 ----------------------------------
const search = (value) => {
  if (!value) return (searchChartsGroupList.value = null)
  // newGroup中滤除不包含value的chart，需要注意的是不能直接修改原数组，否则会影响到原数组
  const newGroups = []
  for (const group of chartsGroupList.value) {
    if (group.charts.some((chart) => chart.name.toLowerCase().includes(value.toLowerCase()))) {
      newGroups.push({ ...group, charts: [...group.charts] })
    }
  }
  newGroups.forEach((group) => {
    group.charts = [...group.charts]
    group.charts = group.charts.filter((chart) => {
      return chart.name.toLowerCase().includes(value.toLowerCase())
    })
  })

  searchChartsGroupList.value = newGroups
}

// ---------------------------------- 依赖注入 ----------------------------------
provide('media', props.media)
provide('defaultColor', props.defaultColor)
provide('getColor', props.getColor)
provide('$on', props.subscribe)
provide('$off', props.unsubscribe)
provide('experimentPrefix', props.experimentPrefix)

console.log(props.experimentPrefix)
</script>

<style lang="scss" scoped>
.charts-pannel {
  @apply sticky px-4 py-3 bg-higher flex justify-between top-0 border-b;
  z-index: 30;
}
</style>
