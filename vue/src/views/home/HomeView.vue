<template>
  <div class="p-6 flex flex-col gap-5 text-dimmer border-b">
    <h1 class="text-2xl font-semibold text-default">{{ projectStore.name }}</h1>
    <!-- <p>{{ projectStore.description }}</p> -->
    <!-- 项目创建时间、最近运行的时间、总实验数量 -->
    <div class="w-80 flex flex-col gap-4">
      <div class="flex justify-between">
        <div class="grow">{{ $t('home.overview.create') }}</div>
        <p class="w-36 whitespace-nowrap">{{ createTime }}</p>
      </div>
      <div class="flex justify-between">
        <div class="grow">{{ $t('home.overview.latest') }}</div>
        <p class="w-36 whitespace-nowrap">{{ updateTime }}</p>
      </div>
      <div class="flex justify-between">
        <div class="grow">{{ $t('home.overview.total') }}</div>
        <p class="w-36 whitespace-nowrap">{{ projectStore.sum }}</p>
      </div>
    </div>
  </div>
  <div class="p-6">
    <h2 class="text-xl font-semibold mb-4">{{ $t('home.list.title') }}</h2>
    <!-- 实验表格 -->
    <SLTable :column="column" :data="experiments_table" high-light v-if="tags">
      <template v-slot:name="{ row }">
        <ExperimentName :name="row.name" :id="row.experiment_id" :color="row.color" />
      </template>
      <template v-slot:status="{ row }">
        <SLStatusLabel :id="row.experiment_id" :status="row.status" />
      </template>
      <template v-slot:create="{ row }">
        {{ transTime(convertUtcToLocal(row.create_time)) }}
      </template>
      <template v-for="item in configs" :key="item.key" v-slot:[item.key]="{ row }">
        {{ row.config[item.key] || '-' }}
      </template>
    </SLTable>
  </div>
</template>

<script setup>
/**
 * @description: 首页视图，列出项目的基本信息
 * @file: HomeView.vue
 * @since: 2023-12-04 19:36:21
 **/
import { useProjectStore } from '@swanlab-vue/store'
import { formatTime } from '@swanlab-vue/utils/time'
import { computed, ref } from 'vue'
import SLStatusLabel from '@swanlab-vue/components/SLStatusLabel.vue'
import ExperimentName from './components/ExperimentName.vue'
import { transTime, convertUtcToLocal } from '@swanlab-vue/utils/time'
import SLTable from '@swanlab-vue/components/SLTable.vue'
import { t } from '@swanlab-vue/i18n'
import http from '@swanlab-vue/api/http'

const projectStore = useProjectStore()

// ---------------------------------- 在此处处理项目创建时间、运行时间和总实验数量 ----------------------------------
const createTime = computed(() => formatTime(projectStore.createTime))
const updateTime = computed(() => formatTime(projectStore.updateTime))

// ---------------------------------- 表格配置 ----------------------------------

const column = ref([
  {
    title: t('home.list.table.header.name'),
    slot: 'name',
    style: 'px-4',
    width: '200'
    // fixed: true
  },
  {
    title: t('home.list.table.header.status'),
    slot: 'status'
  },
  {
    title: t('home.list.table.header.create'),
    slot: 'create'
  }
])

// 遍历 projectStore.experiments，添加配置
const configs = []
;(() => {
  // 寻找需要增加的表头
  projectStore.experiments.map((item) => {
    Object.entries(item.config).forEach(([key]) => {
      // 如果这个key已经存在configs中，跳过
      if (configs.some((config) => config.title === key)) {
        return
      }
      configs.push({
        key,
        slot: key,
        title: key
      })
    })
  })
  column.value.push(...configs)
})()

// ---------------------------------- 表格数据，同时还有tag的表头处理 ----------------------------------

// 表格体数据
const experiments_table = computed(() => {
  return projectStore.experiments.map((expr) => {
    const summary = summaries.value[expr.name]
    if (!summary) return {}
    Promise.all(
      Object.keys(summary).map(async (key) => {
        expr[await hashString(key)] = summary[key]
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
        projectStore.experiments.forEach((experiment) => {
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
    column.value.push(...tags.value)
    // 保存tag总结数据
    summaries.value = data.summaries
  })

// 哈希处理 key 避免和关键字重复
async function hashString(inputString) {
  const encoder = new TextEncoder()
  const data = encoder.encode(inputString)

  const buffer = await crypto.subtle.digest('SHA-256', data)
  const hashArray = Array.from(new Uint8Array(buffer))
  const hashHex = hashArray.map((byte) => byte.toString(16).padStart(2, '0')).join('')

  return hashHex
}
</script>

<style lang="scss" scoped>
.experiments-table {
  @apply border;
  tr {
    &:first-child {
      @apply bg-dimmer;
    }
    &:not(:first-child) {
      @apply border-t;
    }
  }

  th {
    @apply text-left;
  }

  th,
  td {
    @apply px-5 py-2.5;

    &:not(:last-child) {
      @apply border-r;
    }
  }
}
</style>
