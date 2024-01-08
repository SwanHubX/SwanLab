<template>
  <div class="w-full h-full px-7 py-6 relative overflow-hidden">
    <SLSearch class="max-w-[400px]" @input="search" dealy="100" />
    <section class="log-container">
      <div class="log-area" ref="logAreaRef" v-if="logs">
        <!-- 运行日志 -->
        <div class="log-line" v-for="line in filteredLogs" :key="line">
          <!-- 行号 -->
          <span class="w-8 text-right flex-shrink-0 text-dimmest select-none">{{ line.index }}</span>
          <!-- 日志内容 -->
          <!-- 如果没有搜索内容/不含搜索内容 -->
          <span v-show="!searchValue || !line.value.includes(searchValue)">{{ line.value }} </span>
          <!-- 如果有搜索内容且含有搜索内容 -->
          <p class="flex" v-show="searchValue && line.value.includes(searchValue)">
            <span>{{ line.value.substring(0, line.value.indexOf(searchValue)) }}</span>
            <!-- 高亮展示 -->
            <span class="bg-negative-default text-white-default">{{ searchValue }}</span>
            <span>{{ line.value.substring(line.value.indexOf(searchValue) + searchValue.length) }}</span>
          </p>
        </div>
        <!-- 错误日志 -->
        <div class="log-line text-negative-default" v-for="(line, index) in errorLogs" :key="line">
          <!-- 行数 -->
          <span>{{ logs.length + index }}</span>
          <!-- 日志内容 -->
          <span>{{ line }}</span>
        </div>
      </div>
      <div class="flex h-full items-center justify-center" v-else>
        <SLLoding />
      </div>
    </section>
    <SLOperation class="absolute top-10 right-10">
      <span>...</span>
      <template #main> nihao </template>
    </SLOperation>
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
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import SLLoding from '@swanlab-vue/components/SLLoading.vue'
import SLOperation from '@swanlab-vue/components/SLOperation.vue'
import SLSearch from '@swanlab-vue/components/SLSearch.vue'
import { computed } from 'vue'

const logAreaRef = ref()

// ---------------------------------- 系统相关 ----------------------------------

const experimentStore = useExperimentStroe()
const id = experimentStore.id

// ---------------------------------- 获取日志和日志相关数据 ----------------------------------

// 所有日志
const logs = ref()

// 分开行号和内容之后的日志
const filteredLogs = computed(() => {
  return logs.value.map((line) => {
    return {
      index: line.substring(0, line.indexOf(' ')),
      value: line.substring(line.indexOf(' '))
    }
  })
})

// 错误日志
const errorLogs = ref([])
;(async function () {
  // 获取日志
  const { data } = await http.get(`/experiment/${id}/recent_log`)
  // 设置日志
  logs.value = data.logs
  if (data.error) errorLogs.value = data.error
  addTaskToBrowserMainThread(() => {
    // 滚动到底部
    logAreaRef.value.scrollTop = logAreaRef.value.scrollHeight
  })
})()

// ---------------------------------- 搜索 ----------------------------------

const searchValue = ref('')

const search = (value) => {
  searchValue.value = value
}
</script>

<style lang="scss" scoped>
.log-container {
  @apply bg-higher w-full h-full rounded p-4;
  font-size: 13px;
  line-height: 16px;
  font-family: 'JetBrains Mono', monospace;
  letter-spacing: 0.1px;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  .log-area {
    @apply overflow-auto h-full;
    &::-webkit-scrollbar-track {
      background: transparent;
    }
  }

  .log-line {
    @apply flex gap-2 whitespace-pre-wrap;
    span {
      @apply block;
    }
  }
}
</style>
