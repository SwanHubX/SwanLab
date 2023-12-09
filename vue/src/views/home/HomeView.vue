<template>
  <div class="p-6 flex flex-col gap-5 text-dimmer border-b">
    <h1 class="text-2xl font-semibold text-default">{{ projectStore.name }}</h1>
    <p>{{ projectStore.description }}</p>
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
    <table class="experiments-table">
      <tr>
        <th class="w-52">{{ $t('home.list.table.header.name') }}</th>
        <th class="w-36">{{ $t('home.list.table.header.status') }}</th>
        <th class="w-36">{{ $t('home.list.table.header.create') }}</th>
      </tr>
      <tr v-for="experiment in projectStore.experiments" :key="experiment.experiment_id">
        <!-- 实验名称 -->
        <td>
          <ExperimentName :name="experiment.name" :id="experiment.experiment_id" :color="experiment.color" />
        </td>
        <!-- 实验状态 -->
        <td>
          <StatusLabel :name="experiment.name" :id="experiment.experiment_id" :status="experiment.status" />
        </td>
        <!-- 创建时间 -->
        <td>{{ transTime(convertUtcToLocal(experiment.create_time)) }}</td>
      </tr>
    </table>
  </div>
</template>

<script setup>
/**
 * @description: 首页视图，列出项目的基本信息
 * @file: HomeView.vue
 * @since: 2023-12-04 19:36:21
 **/
import { useProjectStore } from '@swanlab-vue/store'
import { formatTime } from '@swanlab-vue/utils/common'
import { computed } from 'vue'
import StatusLabel from '@swanlab-vue/components/StatusLabel.vue'
import ExperimentName from './components/ExperimentName.vue'
import { transTime, convertUtcToLocal } from '@swanlab-vue/utils/time'

const projectStore = useProjectStore()

// ---------------------------------- 在此处处理项目创建时间、运行时间和总实验数量 ----------------------------------
const createTime = computed(() => formatTime(projectStore.createTime))
const updateTime = computed(() => formatTime(projectStore.updateTime))
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
