<template>
  <section class="chart-container" @mouseenter="handleMouseEnter" @mouseleave="handleMouseLeave">
    <template v-if="status === 'success'">
      <!-- 图表标题 -->
      <p class="text-center font-semibold">{{ title(chart.source[0]) }}</p>
      <!-- 图表相关控制按钮 -->
      <div class="chart-pannel" v-if="hover"></div>
      <component />
    </template>
    <!-- 错误 -->
    <div class="flex flex-col justify-center grow text-dimmer gap-2" v-else-if="status === 'error'">
      <SLIcon class="mx-auto h-5 w-5" icon="error" />
      <p class="text-center text-xs">此图表无法被正确显示</p>
    </div>
    <!-- 加载中 -->
    <div class="h-full flex items-center justify-center" v-else>
      <SLLoading />
    </div>
  </section>
</template>

<script setup>
/**
 * @description: 图表容器，用于包裹并渲染图表组件
 * @file: ChartContainer.vue
 * @since: 2023-12-12 21:01:26
 **/
import { inject } from 'vue'
import { onUnmounted } from 'vue'
import { ref } from 'vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import SLLoading from '@swanlab-vue/components/SLLoading.vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
const props = defineProps({
  chart: {
    type: Object,
    required: true
  }
})

const cid = props.chart._cid
const source = props.chart.source

// ---------------------------------- 判断是否处于hover状态 ----------------------------------
const hover = ref(false)
const handleMouseEnter = () => (hover.value = true)
const handleMouseLeave = () => (hover.value = false)

// ---------------------------------- 控制chart显示状态 ----------------------------------
const status = ref('loading')

const chartComponent = (type) => {}

// ---------------------------------- 订阅 ----------------------------------
const $off = inject('$off')
inject('$on')(source, cid, (tag, data, error) => {
  return new Promise((resolve, reject) => {
    if (error) {
      status.value = 'error'
      reject(error)
    } else {
      status.value = 'success'
      // console.log(tag, data)
      // 完成渲染
      addTaskToBrowserMainThread(() => {
        console.log('渲染图表: ', props.chart)
      })
      resolve(data)
    }
  })
})

// 卸载时取消订阅
onUnmounted(() => {
  $off(source, cid)
})

// ---------------------------------- 图表标题控制 ----------------------------------
const title = (t) => {
  // TODO 多数据时可能有不同的显示方式
  return t
}
</script>

<style lang="scss" scoped>
.chart-container {
  @apply w-full h-72 border rounded relative overflow-hidden;
  @apply px-3 py-4;
  @apply flex-col flex justify-between;
}

.chart-pannel {
  @apply absolute top-1 right-2   h-4;
  @apply flex justify-end gap-2;
}
</style>
