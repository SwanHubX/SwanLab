<template>
  <div class="w-full h-full px-7 py-6">
    <section class="log-container">
      <div class="log-area" ref="logAreaRef" v-if="logs">
        <div class="log-line" v-for="line in logs" :key="line">
          <!-- 行号 -->
          <span>{{ line.substring(0, line.indexOf(' ')) }}</span>
          <!-- 日志内容 -->
          <span>{{ line.substring(line.indexOf(' ')) }}</span>
        </div>
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
const logAreaRef = ref()
// ---------------------------------- 系统相关 ----------------------------------

const experimentStore = useExperimentStroe()
const id = experimentStore.id

// ---------------------------------- 获取日志和日志相关数据 ----------------------------------
// 所有日志
const logs = ref()
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
    @apply flex gap-2;
    span {
      @apply block;
    }
    span:first-child {
      @apply w-8 text-right flex-shrink-0 text-dimmest select-none;
    }
  }
}
</style>
