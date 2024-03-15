<template>
  <HomeLayout>
    <template #top>
      <!-- 项目创建时间、最近运行的时间、总实验数量 -->
      <div class="w-80 flex flex-col gap-4 my-5 text-dimmer">
        <div class="flex justify-between">
          <div class="grow">{{ $t('home.overview.create') }}</div>
          <p class="w-36 whitespace-nowrap">{{ createTime }}</p>
        </div>
        <div class="flex justify-between">
          <div class="grow">{{ $t('home.overview.latest') }}</div>
          <p class="w-36 whitespace-nowrap">{{ updateTime }}</p>
        </div>
        <div class="flex justify-between">
          <div class="grow">{{ $t('home.overview.logdir') }}</div>
          <p class="w-36 whitespace-nowrap">{{ projectStore.logdir }}</p>
        </div>
      </div>
    </template>
    <div class="pt-5 flex flex-col flex-grow flex-shrink-0">
      <!-- 实验列表+实验统计 -->
      <div class="flex mx-6 mb-4 items-center gap-1.5 flex-shrink-0">
        <h2 class="text-xl font-semibold">{{ $t('home.list.title') }}</h2>
        <span v-if="total" class="bg-highest text-default rounded-full px-3">{{ total }}</span>
      </div>
      <TableBar
        class="pb-4 pt-1 px-5 flex-shrink-0"
        :table-head="tableHead"
        :table-body="tableBody"
        :searchText="searchText"
        :checked="onlySummary"
        @update:checked="onlySummary = $event"
        @update:searchText="searchText = $event"
      />
      <div class="w-full pb-10 flex flex-col overflow-scroll flex-grow basis-0" v-if="tags && init">
        <!-- 实验表格 -->
        <ExprTable sticky-header class="dashboard-table" :column="tableHead" :data="tableBody" last-row-gradient>
          <template v-slot:name="{ row }">
            <ExperimentName :name="row.name" :id="row.id" :color="getExperimentColor(row)" />
          </template>
          <template v-slot:status="{ row }">
            <SLStatusLabel :name="row.name" :status="row.status" :url="'/experiment/' + row.id">
              {{ $t('experiment.status.' + row.status) }}
            </SLStatusLabel>
          </template>
          <template v-slot:create="{ row }">
            {{ transTime(convertUtcToLocal(row.create_time)) }}
          </template>
          <template v-slot:duration="{ row }">
            {{ getDuration(row) }}
          </template>
          <template v-for="item in configs" :key="item.slot" v-slot:[item.slot]="{ row }">
            {{ row.config[item.slot]?.value != undefined ? row.config[item.slot]?.value : '-' }}
          </template>
        </ExprTable>
      </div>
      <EmptyTable v-else-if="experiments.length === 0 && init" />
      <div v-else>
        <div class="w-full flex" v-for="index in 3" :key="index">
          <div v-for="i in 5" :key="i" class="skeleton w-full py-3 m-2"></div>
        </div>
      </div>
    </div>
  </HomeLayout>
</template>

<script setup>
/**
 * @description: 首页视图，列出项目的基本信息
 * @file: HomeView.vue
 * @since: 2023-12-04 19:36:21
 **/
import { useProjectStore } from '@swanlab-vue/store'
import { formatTime, getDuration } from '@swanlab-vue/utils/time'
import { computed, ref, onMounted } from 'vue'
import SLStatusLabel from '@swanlab-vue/components/SLStatusLabel.vue'
import ExperimentName from './components/ExperimentName.vue'
import { transTime, convertUtcToLocal } from '@swanlab-vue/utils/time'
import { t } from '@swanlab-vue/i18n'
import http from '@swanlab-vue/api/http'
import ExprTable from './components/ExprTable.vue'
import EmptyTable from './components/EmptyTable.vue'
import TableBar from './components/TableBar.vue'
import { formatNumber2SN } from '@swanlab-vue/utils/common'

const projectStore = useProjectStore()
const experiments = computed(() => {
  // 在最前面判断项目信息是否存在，不存在则是后端未开启/没有项目
  if (typeof projectStore.experiments === 'undefined') {
    return []
  }
  return projectStore.experiments
})

// ---------------------------------- 图表延迟 ----------------------------------

const init = ref(false)
onMounted(() => {
  setTimeout(() => {
    init.value = true
  }, 500)
})

// ---------------------------------- 在此处处理项目创建时间、运行时间和总实验数量 ----------------------------------

const createTime = computed(() => formatTime(projectStore.createTime))
const updateTime = computed(() => formatTime(projectStore.updateTime))
const total = computed(() => projectStore.experiments?.length || 0)

// ---------------------------------- 表格配置 ----------------------------------

const column = ref([
  {
    title: t('home.list.table.header.name'),
    slot: 'name',
    style: 'px-6',
    width: 300,
    border: true
  },
  {
    title: t('home.list.table.header.status'),
    slot: 'status',
    type: 'status'
  },
  {
    title: t('home.list.table.header.create'),
    slot: 'create',
    type: 'create_time'
  },
  {
    title: t('home.list.table.header.duration'),
    slot: 'duration',
    type: 'duration'
  }
])

// 遍历 experiments，添加配置
const configs = ref([])
;(() => {
  // 寻找需要增加的表头
  experiments.value.map((item) => {
    Object.entries(item.config).forEach(([key]) => {
      // 如果这个key已经存在configs中，跳过
      if (configs.value.some((config) => config.title === key)) {
        return
      }
      configs.value.push({
        type: 'config',
        slot: key,
        title: key
      })
    })
  })
  // column.value.push(...configs.value)
})()

// ---------------------------------- 表格数据，同时还有tag的表头处理 ----------------------------------

// 表格体数据
const experiments_table = computed(() => {
  return experiments.value.map((expr) => {
    const summary = summaries.value[expr.name]
    if (!summary) return {}
    Promise.all(
      Object.keys(summary).map(async (key) => {
        expr[await hashString(key)] = isNaN(summary[key]) ? summary[key] : formatNumber2SN(summary[key])
      })
    )
    return expr
  })
})

// 项目里面的所有 tag 项，undefined 表示还没有初始化完，这个时候不加载表格
const tags = ref()
const summaries = ref({})
http
  .get('/project/summaries', {
    params: {
      // 传递前端显示的所有实验名称，使用字符串格式，每个实验名称之间使用逗号连接
      experiment_names: (() => {
        let experiment_names = []
        experiments.value.forEach((experiment) => {
          experiment_names.push(experiment.name)
        })
        return experiment_names.join(',')
      })()
    }
  })
  .then(async ({ data }) => {
    // 增加tag对应的表头
    tags.value = await Promise.all(
      data.tags.map(async (tag) => {
        const key = await hashString(tag)
        return { key, title: tag }
      })
    )
    // 保存tag总结数据
    summaries.value = data.summaries
  })

// 哈希处理 key 避免和关键字重复
async function hashString(inputString) {
  // const encoder = new TextEncoder()
  // const data = encoder.encode(inputString)

  // const buffer = await crypto.subtle.digest('SHA-256', data)
  // const hashArray = Array.from(new Uint8Array(buffer))
  // const hashHex = hashArray.map((byte) => byte.toString(16).padStart(2, '0')).join('')

  // return hashHex
  return 'swanlab-overview-table-key' + inputString
}

// ---------------------------------- 颜色选择 ----------------------------------

const getExperimentColor = (experiment) => {
  return experiment.light
}
// ---------------------------------- 筛选 ----------------------------------

// 是否开启了表格筛选
const onlySummary = ref(false)

// 动态修改表头
const tableHead = computed(() => {
  if (!onlySummary.value) {
    return [...column.value, ...configs.value, ...(tags.value || [])]
  }

  return [column.value[0], ...tags.value]
})

// ---------------------------------- 查找 ----------------------------------

const searchText = ref('')

const tableBody = computed(() => {
  if (!searchText.value) {
    return experiments_table.value
  }
  return experiments_table.value.filter((item) => {
    return item.name.toLowerCase().includes(searchText.value.toLowerCase())
  })
})
</script>

<style lang="scss" scoped>
// 最后一行每一个单元格右边框添加伪元素
.dashboard-table {
  @apply border-x-0 border-b-0;
}
.grow {
  @apply text-dimmest;
}

@-webkit-keyframes skeleton-ani {
  0% {
    left: 0;
  }

  to {
    left: 100%;
  }
}

@keyframes skeleton-ani {
  0% {
    left: 0;
  }

  to {
    left: 100%;
  }
}

.skeleton {
  background-image: linear-gradient(-45deg, #f5f5f5 40%, #fff 55%, #f5f5f5 63%);
  list-style: none;
  background-size: 400% 100%;
  background-position: 100% 50%;
  animation: skeleton-animation 2s ease infinite;
}

@keyframes skeleton-animation {
  0% {
    background-position: 100% 50%;
  }

  100% {
    background-position: 0 50%;
  }
}
</style>
