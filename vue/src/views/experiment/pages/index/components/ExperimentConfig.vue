<template>
  <div class="pr-5" icon="config" :title="$t('experiment.index.config.title')" :retract="false">
    <div class="pl-6 w-full grid lg:grid-cols-2 lg:gap-10">
      <div class="pt-4">
        <div class="flex items-center pb-4">
          <p class="font-semibold pr-2">{{ $t('experiment.index.config.detail') }}</p>
          <SLHelp document="https://geektechstudio.feishu.cn/wiki/EFi3wuACGiEWlLki5aDcQiSpngg">{{
            $t('experiment.index.config.help.config')
          }}</SLHelp>
        </div>
        <DataTable table-border :column="column" :data="configs" title="Config" />
      </div>
      <div class="pt-4">
        <div class="flex items-center pb-4">
          <p class="font-semibold pr-2">{{ $t('experiment.index.config.summarize') }}</p>
          <SLHelp document="https://geektechstudio.feishu.cn/wiki/TudNwOSMyihFetky7l5cTI8UnJf"
            >{{ $t('experiment.index.config.help.summary') }}
          </SLHelp>
        </div>
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
import SLHelp from '@swanlab-vue/components/SLHelp.vue'
import http from '@swanlab-vue/api/http'
import { ref } from 'vue'
import { useExperimentStroe } from '@swanlab-vue/store'
import { computed } from 'vue'
import { formatNumber2SN } from '@swanlab-vue/utils/common'

const experiment = ref(useExperimentStroe().experiment)

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
      value: experiment.value.config[key].value
    })
  }
  return configs
})

// ---------------------------------- 获取实验的总结数据 ----------------------------------

const summaries = ref([])
const experimentId = ref(useExperimentStroe().id)
http
  .get(`/experiment/${experimentId.value}/summary`)
  .then(({ data }) => {
    summaries.value = (data.summaries || []).map((item) => {
      return {
        key: item.key,
        value: formatNumber2SN(item.value)
      }
    })
  })
  .catch((error) => {
    console.error(error)
  })
</script>

<style lang="scss" scoped></style>
