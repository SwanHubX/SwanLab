<template>
  <div class="w-full h-full px-7 py-6">
    <section class="log-area">
      <div class="overflow-auto h-full" ref="logAreaRef">
        <div class="log-line" v-for="(line, index) in logs" :key="line">
          <!-- 行数 -->
          <span class="select-none">{{ index + 1 }}</span>
          <!-- 日志内容 -->
          <span>{{ line }}</span>
        </div>
        <div class="log-line text-negative-default" v-for="(line, index) in errorLogs" :key="line">
          <!-- 行数 -->
          <span class="select-none">{{ logs.length + index + 1 }}</span>
          <!-- 日志内容 -->
          <span>{{ line }}</span>
        </div>
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
const logAreaRef = ref()
// ---------------------------------- 系统相关 ----------------------------------

const experimentStore = useExperimentStroe()
const id = experimentStore.id

// ---------------------------------- 获取日志和日志相关数据 ----------------------------------
// 所有日志
const logs = ref([])
// 错误日志
const errorLogs = ref([])
;(async function () {
  // 获取日志
  const { data } = await http.get(`/experiment/${id}/log`)
  // 设置日志
  logs.value = data.logs
  addTaskToBrowserMainThread(() => {
    // 滚动到底部
    logAreaRef.value.scrollTop = logAreaRef.value.scrollHeight
  })
})()
</script>

<style lang="scss" scoped>
.log-area {
  @apply bg-dimmer w-full h-full rounded p-4;
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
