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
    <TextModule
      class="text-table"
      :data="original_data[tag]"
      :texts="texts[tag]"
      :tag="tag"
      @getText="getText"
      v-for="tag in source"
      :key="tag"
    />
    <!-- 放大效果弹窗 -->
    <SLModal class="pb-4 overflow-hidden" max-w="-1" v-model="isZoom">
      <TextModule
        :data="original_data[tag]"
        :texts="texts[tag]"
        :tag="tag"
        @getText="getText"
        v-for="tag in source"
        :key="tag"
        modal
      />
    </SLModal>
  </template>
</template>

<script setup>
/**
 * @description: 文字图标
 * @file: TextChart.vue
 * @since: 2024-02-20 20:04:09
 **/
import { ref, computed } from 'vue'
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

// 文本内容，数据结构：{ tag: Array<array || string> }
const texts = ref({})

/**
 * 通过 log 中记录的文件名获取文本内容
 * @param {string} tag 标签名
 * @param {int} currentPage 数据在 list 中的索引号 => step 对应的索引，指出了数据在数组中的真实位置
 */
const getText = async (tag, currentPage) => {
  // 已经有数据了，不需要重复请求
  if (texts.value[tag][currentPage] != null) return
  // 即训练时的 log 内容，保存了文本文件名称
  const path = original_data.value[tag].list[currentPage].data
  const res = await http.get('/media/text', {
    params: {
      // 如果 path 是数组，则拼接成字符串 => 该接口通过 params 获取数组
      path: Array.isArray(path) ? path.join(',') : path,
      tag,
      run_id: run_id.value
    }
  })
  texts.value[tag][currentPage] = res.data.text
}

// ---------------------------------- 渲染、重渲染功能 ----------------------------------

// 渲染
const render = (data) => {
  // 保存原始数据，其中主要数据结构：{tag: {list: [{data: 'xxx', more: {caption: 'xxx'}}]}}
  original_data.value = data
  // 遍历 source 获取所有 tag 第一页的文本数据
  for (let index in source.value) {
    const tag = source.value[index]
    // texts 是对象，属性名为 tag，属性值为数组，数组长度为每个 tag 所有数据的总数，未获取的数据位置用 null 填充
    texts.value[tag] = Array(data[tag].list.length).fill(null)
    // 第二参数为 0，表示获取 tag 第一页的文本
    getText(tag, 0)
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
