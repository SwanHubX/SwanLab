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
    <section class="log-container" @scroll="handleScroll">
      <!-- 用于撑开容器，保证滚动条高度 -->
      <div class="log-area" :style="{ height: areaHeight + 'px' }" v-if="logs">
        <!-- 日志画板 -->
        <div class="log-board" :style="{ top: computeTop + 'px' }">
          <!-- 运行日志 -->
          <div
            class="log-line"
            v-for="(line, index) in lines.slice(range[0], range[1] + 1)"
            :key="line.value + range[0] + index"
          >
            <!-- 行号 -->
            <span
              class="w-8 text-right flex-shrink-0 text-dimmest select-none"
              :class="{ 'text-negative-default': line.error }"
              :style="{ width: indexWidth }"
              >{{ line.index }}</span
            >
            <!-- 正常日志 -->
            <div v-if="!line.error">
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
            <div class="text-negative-default" v-else>{{ line.value }}</div>
          </div>
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

// ---------------------------------- 系统相关 ----------------------------------

const experimentStore = useExperimentStore()
const id = experimentStore.id

// ---------------------------------- 获取日志和日志相关数据 ----------------------------------

// 所有日志
const logs = ref()

// 分开行号和内容之后的日志
const lines = computed(() => {
  // 正常日志内容
  const data = logs.value?.map((line, index) => {
    const content = line
    const isTarget = searchValue.value !== '' && content.toLowerCase().includes(searchValue.value)
    return {
      isTarget,
      index: index + 1,
      value: isTarget ? splitStringBySearch(content, searchValue.value) : content
    }
  })
  // 正常内容最后一行的行号
  const lastLineIndex = getLastLineIndex(data)
  // 处理错误日志内容
  errorLogs.value.map((line, index) => {
    return data.push({
      index: lastLineIndex + index + 1,
      value: line,
      error: true
    })
  })

  return data
})

// 行号宽度
const indexWidth = computed(() => {
  let index = 1
  const l = lines.value?.length
  if (l) {
    index = lines.value[l - 1].index === '' ? lines.value[l - 2].index : lines.value[l - 1].index
  }
  if (!index) return '0px'
  return 8 * index.length + 'px'
})

// 获取正常内容最后一行的行号
const getLastLineIndex = (data) => {
  const maxIndex = data?.reduce((max, line) => {
    if (!line.index) return max
    return Number(line.index) > Number(max) ? line.index : max
  }, 0)
  return Number(maxIndex)
}

const download = computed(() => {
  return logs.value?.join('\n') + errorLogs.value.join('\n')
})

const splitStringBySearch = (target, substring) => {
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
const lineHeight = ref(16)
// 最大额外渲染行数
const addition = 30
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

/**
 * 计算 log 行高
 * 没有使用计算属性，因为担心用户手动更改缩放等情况，不会触发响应式计算
 * 采用函数，对某些事件进行监听后调用以适应不同情况下的行高
 * 一般行高是 16px，基本不会出问题，所以目前只在初始化时计算，而没做另外的事件监听
 */
const computeLineHeight = () => {
  const e = document.querySelector('.log-line')
  if (!e) return 0
  const line_styles = window.getComputedStyle(e)
  lineHeight.value = parseFloat(line_styles.getPropertyValue('line-height'))
}

const computeTop = computed(() => {
  return range.value[0] * lineHeight.value
})

/**
 * 初始化：
 * 1. 获取日志数据
 * 2. 计算日志行高
 * 3. 计算日志渲染范围
 */
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

  // 保证渲染后再获取行高，再计算渲染范围
  addTaskToBrowserMainThread(() => {
    // 计算
    computeRange(document.querySelector('.log-container'))
    computeLineHeight()

    // 判断是否出现滚动条，如果有滚动条，滚动到底部
    if (document.querySelector('.log-container').scrollHeight > document.querySelector('.log-container').clientHeight) {
      document.querySelector('.log-container').scrollTop = document.querySelector('.log-container').scrollHeight
    }
  })
})
</script>

<style lang="scss" scoped>
.log-container {
  @apply bg-higher w-full rounded p-4 h-[76vh] overflow-auto relative;
  font-size: 13px;
  line-height: 16px;
  font-family: 'JetBrains Mono', monospace;
  letter-spacing: 0.1px;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;

  .log-area {
    @apply relative h-full break-all;
    &::-webkit-scrollbar-track {
      background: transparent;
    }

    .log-board {
      @apply absolute pb-5 w-full;
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
