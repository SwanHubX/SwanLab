<template>
  <!-- 项目信息 -->
  <div class="p-6 flex flex-col gap-5 text-dimmer border-b relative">
    <!-- 删除项目 -->
    <SLDelete class="absolute top-5 right-4" type="project" @confirm="deleteProject" :disabled="hasRunning" />
    <div class="flex gap-3">
      <h1 class="text-2xl font-semibold text-default">{{ projectStore.name }}</h1>
      <ConfigEditor type="project" @modify="modifyProject" />
    </div>
    <p v-if="projectStore.description">{{ projectStore.description }}</p>
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
    <SLTable :column="column" :data="experiments_table" v-if="tags">
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
    <EmptyTable v-else-if="experiments.length === 0" />
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
import { computed, ref, inject } from 'vue'
import SLStatusLabel from '@swanlab-vue/components/SLStatusLabel.vue'
import ExperimentName from './components/ExperimentName.vue'
import { transTime, convertUtcToLocal } from '@swanlab-vue/utils/time'
import { t } from '@swanlab-vue/i18n'
import http from '@swanlab-vue/api/http'
import ConfigEditor from '@swanlab-vue/components/config-editor/ConfigEditor.vue'
import SLTable from '@swanlab-vue/components/table'
import SLDelete from '@swanlab-vue/components/SLDelete.vue'
import { message } from '@swanlab-vue/components/message'
import EmptyTable from './components/EmptyTable.vue'

const projectStore = useProjectStore()
const experiments = computed(() => {
  // 在最前面判断项目信息是否存在，不存在则是后端未开启/没有项目
  if (typeof projectStore.experiments === 'undefined') {
    show_error(3500)
    return []
  }
  return projectStore.experiments
})

const show_error = inject('show_error')

// ---------------------------------- 在此处处理项目创建时间、运行时间和总实验数量 ----------------------------------

const createTime = computed(() => formatTime(projectStore.createTime))
const updateTime = computed(() => formatTime(projectStore.updateTime))

// ---------------------------------- 表格配置 ----------------------------------

const column = ref([
  {
    title: t('home.list.table.header.name'),
    slot: 'name',
    style: 'px-4',
    width: 300,
    border: true
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

// 遍历 experiments，添加配置
const configs = []
;(() => {
  // 寻找需要增加的表头
  experiments.value.map((item) => {
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
  return experiments.value.map((expr) => {
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
    column.value.push(...tags.value)
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

// ---------------------------------- 修改项目信息 ----------------------------------

const modifyProject = async (newV, hideModal) => {
  const { data } = await http.patch('/project/update', newV)
  projectStore.setProject(data.project)
  hideModal()
}

// ---------------------------------- 删除项目 ----------------------------------

/**
 * 是否可以删除项目
 * 如果有实验正在进行中，直接不显示删除按钮
 */
const hasRunning = computed(() => {
  for (const experiment of experiments.value) {
    if (experiment.status === 0) return true
  }
  return false
})

/**
 * 删除项目
 *
 * success：显示成功，清空状态，并显示空项目错误页
 *
 * fail：在有实验正在进行的时候删除项目
 */
const deleteProject = () => {
  http
    .delete('/project/delete')
    .then(() => {
      message.success('Delete Successfully')
      projectStore.clearProject()
      show_error(3500)
    })
    .catch(({ data }) => {
      message.error('Error deleting project')
      show_error(data.code, data.message)
    })
}
</script>

<style lang="scss" scoped>
.experiments-table {
  @apply border;
  tr {
    &:first-child {
      @apply bg-higher;
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
