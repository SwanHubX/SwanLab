<template>
  <div class="pr-5" icon="config" :title="$t('experiment.index.config.title')" :retract="false">
    <div class="pl-6 w-full grid lg:grid-cols-2 lg:gap-10">
      <div class="pt-4 w-full">
        <DataTable class="w-full" table-border :column="column" :data="configs" title="Config" />
      </div>
      <div class="pt-4">
        <DataTable table-border :column="column" :data="summaries" title="Summary" />
      </div>
    </div>
    <div class="w-full h-6"></div>
  </div>
</template>

<script setup>
/**
 * @description: 实验的配置信息
 * @file: ExperimentConfig.vue
 * @since: 2023-12-11 17:07:31
 **/

import DataTable from './DataTable.vue'
import http from '@swanlab-vue/api/http'
import { ref } from 'vue'
import { useExperimentStore } from '@swanlab-vue/store'
import { computed } from 'vue'
import { formatNumber2SN } from '@swanlab-vue/utils/common'

const experiment = ref(useExperimentStore().experiment)

// 通用表头
const column = [
  {
    title: 'key',
    key: 'key'
  },
  {
    title: 'value',
    key: 'value'
  }
]

// ---------------------------------- 转化实验配置为表格数据 ----------------------------------

const configs = computed(() => {
  const configs = []
  for (const key in experiment.value.config) {
    configs.push({
      key,
      value: experiment.value.config[key].value,
      sort: experiment.value.config[key].sort
    })
  }
  // 按照 sort 排序，然后删除 sort 属性
  return configs
    .sort((a, b) => {
      return a.sort - b.sort
    })
    .map((item) => {
      delete item.sort
      return item
    })
})

// ---------------------------------- 获取实验的总结数据 ----------------------------------

const summaries = ref([])
const experimentId = ref(useExperimentStore().id)
http
  .get(`/experiment/${experimentId.value}/summary`)
  .then(({ data }) => {
    summaries.value = (data.summaries || []).map((item) => {
      return {
        key: item.key,
        value: isNaN(item.value) ? item.value : formatNumber2SN(item.value)
      }
    })
  })
  .catch((error) => {
    console.error(error)
  })
</script>

<style lang="scss" scoped></style>
