<template>
  <div class="w-full h-full px-7 py-6 relative overflow-hidden">
    <FuncBar
      class="pb-6"
      @input="search"
      :content="download"
      :placeholder="$t('experiment.func-bar.placeholder.log')"
      :filename="filename"
      v-if="logs"
    />
    <section class="log-container">
      <div class="log-area" ref="logAreaRef" @scroll="handleScroll" v-if="logs">
        <!-- 运行日志 -->
        <div class="log-line" v-for="line in lines" :key="line">
          <!-- 行号 -->
          <span class="w-8 text-right flex-shrink-0 text-dimmest select-none">{{ line.index }}</span>
          <!-- 日志内容 -->
          <!-- 如果没有搜索内容/不含搜索内容 -->
          <span v-show="!line.isTarget">{{ line.value }} </span>
          <!-- 如果有搜索内容且含有搜索内容 -->
          <span v-show="line.isTarget">
            <span
              v-for="substring in line.value"
              :key="substring"
              :class="{ 'bg-warning-dimmest': substring.toLowerCase() === searchValue }"
            >
              {{ substring }}
            </span>
          </span>
        </div>
        <!-- 错误日志 -->
        <div class="log-line text-negative-default" v-for="(line, index) in errorLogs" :key="line">
          <!-- 行数 -->
          <span class="w-8 text-right flex-shrink-0 select-none">{{ lastLineIndex + index + 1 }}</span>
          <!-- 日志内容 -->
          <span>{{ line }}</span>
        </div>
        <!-- 没有日志 -->
        <div
          class="w-full py-10 flex flex-col items-center text-lg font-semibold"
          v-if="logs.length === 0 && errorLogs.length === 0"
        >
          <SLIcon class="magnifier" icon="search"></SLIcon>
          <span>No Logs</span>
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

import { ref, onMounted } from 'vue'
import http from '@swanlab-vue/api/http'
import { useExperimentStore } from '@swanlab-vue/store'
import SLLoding from '@swanlab-vue/components/SLLoading.vue'
import FuncBar from '../../components/FuncBar.vue'
import { computed } from 'vue'
import { debounce } from '@swanlab-vue/utils/common'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'

const logAreaRef = ref()

// ---------------------------------- 系统相关 ----------------------------------

const experimentStore = useExperimentStore()
const id = experimentStore.id

// ---------------------------------- 获取日志和日志相关数据 ----------------------------------

// 所有日志
const logs = ref()

// 分开行号和内容之后的日志
const lines = computed(() => {
  return logs.value?.map((line) => {
    const noIndex = isNaN(line.substring(0, line.indexOf(' ')))
    const index = noIndex ? null : line.substring(0, line.indexOf(' '))
    const content = noIndex ? line : line.substring(line.indexOf(' ')).trimStart()
    const isTarget = searchValue.value !== '' && content.toLowerCase().includes(searchValue.value)
    return {
      isTarget,
      index,
      value: isTarget ? splitStringBySearch(content, searchValue.value) : content
    }
  })
})

// 正常内容最后一行的行号
const lastLineIndex = computed(() => {
  const maxIndex = lines.value?.reduce((max, line) => {
    if (!line.index) return max
    return Number(line.index) > Number(max) ? line.index : max
  }, 0)
  return Number(maxIndex)
})

const download = computed(() => {
  const log_list = logs.value?.map((log) => {
    const noIndex = isNaN(log.substring(0, log.indexOf(' ')))
    return noIndex ? log : log.substring(log.indexOf(' ')).trimStart()
  })

  return log_list?.join('\n') + errorLogs.value.join('\n')
})

function splitStringBySearch(target, substring) {
  // 使用正则表达式进行大小写不敏感的分割，并保留分割点
  let resultArray = target.split(new RegExp(`(${substring})`, 'i'))

  // 去除数组中的空字符串
  resultArray = resultArray.filter((item) => item !== '')

  // 返回包含分割点的数组
  return resultArray
}

// 错误日志
const errorLogs = ref([])

// ---------------------------------- 搜索 ----------------------------------

const searchValue = ref('')

const search = (value) => {
  searchValue.value = value.toLowerCase()
}

// ---------------------------------- 下载 ----------------------------------

const filename = 'print.log'

// ---------------------------------- 动态渲染 ----------------------------------

// log 渲染范围
const range = ref([0, 0])
// 行高
const lineHeight = ref(0)
// 最大额外渲染行数
const addition = 10
// log 区高度
const areaHeight = computed(() => {
  return lines.value?.length * lineHeight.value
})

const handleScroll = debounce((event) => {
  const e = event.target
  computeRange(e)
}, 100)

// 计算 log 的渲染范围
const computeRange = (e) => {
  const line_height = lineHeight.value
  // 计算应该渲染多少条数据
  const pageSize = Math.ceil(e.clientHeight / line_height)
  // 计算第一条 log 的索引
  const startIndex = Math.floor(e.scrollTop / line_height)
  // 计算最后一条 log 的索引
  const endIndex = startIndex + pageSize

  // 如果距离顶部很近，把顶部到第一条中的log也渲染了
  if (e.scrollTop <= addition * line_height) range.value[0] = 0
  // 如果距离顶部比较远，仅渲染第一条上的 addition 条 log
  else range.value[0] = startIndex - addition

  // 如果距离底部很近，把最后一条到底部的log也渲染了
  if (e.scrollTop + e.clientHeight >= e.scrollHeight - addition * line_height) range.value[1] = lines.value.length - 1
  // 如果距离底部较远，仅渲染最后一条下的 addition 条 log
  else range.value[1] = endIndex + addition
}

// 计算 log 行高
const computeLineHeight = () => {
  const line_styles = window.getComputedStyle(document.querySelector('.log-line'))
  lineHeight.value = parseFloat(line_styles.getPropertyValue('line-height'))
}

onMounted(async () => {
  // 获取日志数据
  const { data } = await http.get(`/experiment/${id}/recent_log`).catch((error) => {
    if (error.data.code === 3404) {
      logs.value = []
    } else {
      console.error(error)
    }
    return
  })
  // 设置日志
  logs.value = data.logs || []
  if (data.error) errorLogs.value = data.error
  // addTaskToBrowserMainThread(() => {
  //   // 滚动到底部
  //   logAreaRef.value.scrollTop = logAreaRef.value.scrollHeight
  // })

  // 保证渲染后再获取行高
  addTaskToBrowserMainThread(computeLineHeight)
})
</script>

<style lang="scss" scoped>
.log-container {
  @apply bg-higher w-full rounded p-4 h-[78vh];
  font-size: 13px;
  line-height: 16px;
  font-family: 'JetBrains Mono', monospace;
  letter-spacing: 0.1px;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  .log-area {
    @apply overflow-auto h-full break-all;
    &::-webkit-scrollbar-track {
      background: transparent;
    }
  }

  .log-line {
    @apply flex gap-2 whitespace-pre-wrap;
  }
}

$duration: 3s;

.magnifier {
  @apply w-10 h-10;
  animation: animloader $duration infinite;
}

@keyframes animloader {
  0% {
    transform: translate(-5px, -5px);
  }
  25% {
    transform: translate(-5px, 5px);
  }
  50% {
    transform: translate(5px, 5px);
  }
  75% {
    transform: translate(5px, -5px);
  }
  100% {
    transform: translate(-5px, -5px);
  }
}
</style>
