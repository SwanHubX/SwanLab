<template>
  <!-- 图表标题 -->
  <p class="text-center font-semibold">{{ title }}</p>
  <!-- 如果图表数据错误 -->
  <div class="flex flex-col justify-center grow text-dimmer gap-2" v-if="error">
    <SLIcon class="mx-auto h-5 w-5" icon="error" />
    <p class="text-center text-xs">
      <!-- 在此处显示错误信息 -->
      {{ error }}
    </p>
  </div>
  <!-- 如果图表数据正确 -->
  <template v-else>
    <!-- 在此处完成图表主体定义 -->
    <TextModule class="text-table" :data="original_data" :texts="texts" :source="source" @getText="getText" />
    <!-- 放大效果弹窗 -->
    <SLModal class="pb-10 overflow-hidden" max-w="-1" v-model="isZoom">
      <TextModule :data="original_data" :texts="texts" :source="source" @getText="getText" modal />
    </SLModal>
  </template>
</template>

<script setup>
/**
 * @description: 文字图标
 * @file: TextChart.vue
 * @since: 2024-02-20 20:04:09
 **/
import { ref, computed, onMounted } from 'vue'
import SLModal from '@swanlab-vue/components/SLModal.vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import TextModule from '../modules/TextModule.vue'
import http from '@swanlab-vue/api/http'
import { useExperimentStore } from '@swanlab-vue/store'

// ---------------------------------- 配置 ----------------------------------

const props = defineProps({
  title: {
    type: String,
    required: true
  },
  chart: {
    type: Object,
    required: true
  },
  index: {
    type: Number,
    required: true
  }
})

const original_data = ref()
const source = ref(props.chart.source)

const experimentStore = useExperimentStore()
const run_id = computed(() => experimentStore.experiment.run_id)

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------

const error = computed(() => {
  if (!props.chart.error) {
    return props.chart.error
  } else if (!original_data.value) {
    return 'No Data'
  }
  return false
})

// ---------------------------------- 获取文本 ----------------------------------

const texts = ref({})

const getText = async (tag, currentPage) => {
  const path = original_data.value[tag].list[currentPage - 1].data
  const res = await http.get('/media/text', {
    params: {
      path: Array.isArray(path) ? path.join(',') : path,
      tag,
      run_id: run_id.value
    }
  })
  texts.value[tag][currentPage - 1] = res.data.text
}

// ---------------------------------- 渲染、重渲染功能 ----------------------------------

// 渲染
const render = (data) => {
  original_data.value = data
  for (let index in source.value) {
    const tag = source.value[index]
    texts.value[tag] = Array(data[tag].list.length).fill(null)
    getText(tag, 1)
  }
}

// 重渲染
const change = (data) => {
  // 将发生更新的 tag 数据保存到原始数据中
  for (let key in data) {
    // original_data.value[key] = data[key] => 这行代码触发不了 props 的响应式，而下面这行可以
    original_data.value[key] = { ...data[key] }
  }
}

// ---------------------------------- 放大功能 ----------------------------------

// 是否放大
const isZoom = ref(false)
// 放大数据
const zoom = () => {
  isZoom.value = true
}

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  render,
  change,
  zoom
})
</script>

<style lang="scss" scoped>
.text-table {
  @apply pt-2 -mx-3 w-[calc(100%+1.5rem)];
}
</style>
