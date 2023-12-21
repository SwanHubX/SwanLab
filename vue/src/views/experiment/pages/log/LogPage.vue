<template>
  <div class="w-full h-full px-7 py-6">
    <section class="log-area">
      <div class="log-line" v-for="(line, index) in logs" :key="line">
        <!-- 行数 -->
        <span class="select-none">{{ index + 1 }}</span>
        <!-- 日志内容 -->
        <span>{{ line }}</span>
      </div>
    </section>
  </div>
</template>

<script setup>
/**
 * @description: 实验-日志页
 * @file: LogPage.vue
 * @since: 2023-12-09 20:40:44
 **/

import { ref } from 'vue'
import http from '@swanlab-vue/api/http'
import { useExperimentStroe } from '@swanlab-vue/store'

// ---------------------------------- 系统相关 ----------------------------------

const experimentStore = useExperimentStroe()
const experiment_id = ref(experimentStore.id)

// ---------------------------------- 获取日志和日志相关数据 ----------------------------------

const logs = ref([])
const total = ref(0)
const current = ref(1)

/**
 * 获取分页对应的日志
 * @param {number} page 分页的页码数
 */
const getLog = async (page) => {
  const res = await http.get(`/experiment/${experiment_id.value}/log`, {
    params: {
      page
    }
  })
  if (res.code != 0) {
    // 错误处理
  }
  total.value = res.data.total
  // 向列表头部添加
  logs.value = [...res.data.logs, ...logs.value]
  console.log(logs.value)
}

// ---------------------------------- 初始化 ----------------------------------

const mini_lines = 2000

const init = async (page) => {
  await getLog(page)
  // 如果行数少于设定值，并且分页未加载完，继续加载下一页
  if (logs.value.length < mini_lines && page < total.value) {
    init(page + 1)
  } else {
    current.value = page
    console.log(current.value)
  }
}

init(current.value)
</script>

<style lang="scss" scoped>
.log-area {
  @apply bg-dimmer w-full h-full overflow-auto rounded p-4;
  font-size: 13px;
  line-height: 16px;
  font-family: 'JetBrains Mono', monospace;
  letter-spacing: 0.1px;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  .log-line {
    @apply flex gap-2;
    span {
      @apply block;
    }
    span:first-child {
      @apply w-8 text-right text-gray-400;
    }
  }
}
</style>
