<template>
  <!-- 图表标题 -->
  <p class="text-center font-semibold">{{ title }}</p>
  <!-- 如果图表数据错误 -->
  <div class="flex flex-col justify-center grow text-dimmer gap-2" v-if="error">
    <SLIcon class="mx-auto h-5 w-5" icon="error" />
    <p class="text-center text-xs">
      <!-- 在此处显示错误信息 -->
      {{ $t('common.chart.charts.image.error', { type: error['data_class'], tag: source[0] }) }}
    </p>
  </div>
  <!-- 如果图表数据正确 -->
  <template v-else>
    <!-- 在此处完成图表主体定义 -->

    <!-- 放大效果弹窗 -->
    <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom"> </SLModal>
  </template>
</template>

<script setup>
/**
 * @description: 图像图表组件
 * @file: ImageChart.vue
 * @since: 2024-02-03 13:24:57
 **/
import SLModal from '@swanlab-vue/components/SLModal.vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import { ref, inject } from 'vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import * as UTILS from './utils'

// ---------------------------------- 配置 ----------------------------------
const props = defineProps({
  title: {
    type: String,
    required: true
  },
  chart: {
    type: Object,
    required: true
  }
})

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------

const error = ref(props.chart.error)

// ---------------------------------- 图表颜色配置 ----------------------------------
// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const colors = inject('colors')
if (!colors) throw new Error('colors is not defined, please provide colors in parent component')
// ---------------------------------- 组件渲染逻辑 ----------------------------------

// ---------------------------------- 数据格式化 ----------------------------------
// 前端映射数据，不包含图像，数据格式：{step: {tag: [{filename: string, caption: string}]}}
const stepsData = {}
// 图像数据缓存, 数据格式：{String<filename>: ImageData}
const imagesData = {}

/**
 * 将后端传回的数据转换为图像数据，存储到imageStepsData中
 * 后端返回的数据格式为：{tag: {..., list:[{data: string or list<strinf>, more: object or list<object>, index: string<number>, create_time: string}]}}
 * @param { Object } data 后端返回的数据
 */
const changeData2Image = (data) => {
  // 遍历数据，将数据转换为图像数据
  for (const tag in data) {
    // 遍历tag下的数据
    for (const item of data[tag].list) {
      if (!stepsData[item.index]) stepsData[item.index] = {}
      // 对当前tag下的数据进行处理,向stepsData中添加数据
      if (!stepsData[item.index][tag]) stepsData[item.index][tag] = []
      // 添加数据,如果data是字符串，则直接添加，如果是数组，则遍历添加
      if (typeof item.data === 'string') {
        stepsData[item.index][tag].push({ filename: item.data, caption: item.more?.caption })
      } else {
        for (let i = 0; i < item.data.length; i++) {
          stepsData[item.index][tag].push({ filename: item.data[i], caption: item.more[i]?.caption })
        }
      }
    }
  }
}

// ---------------------------------- 渲染、重渲染功能 ----------------------------------

// 渲染
const render = (data) => {
  // 数据格式化
  changeData2Image(data)
  console.log('图像数据：', stepsData)
}
// 重渲染
const change = (data) => {
  // 数据格式化
  changeData2Image(data)
}

// ---------------------------------- 放大功能 ----------------------------------
// 是否放大
const isZoom = ref(false)
// 放大数据
const zoom = (data) => {
  isZoom.value = true
  // 放大后图表的高度
  const height = window.innerHeight * 0.6
  addTaskToBrowserMainThread(() => {})
}

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  render,
  change,
  zoom
})
</script>

<style lang="scss" scoped></style>
