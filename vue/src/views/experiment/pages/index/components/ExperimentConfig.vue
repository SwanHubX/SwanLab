<template>
  <ExtendBlock class="pt-2 pr-10" icon="experiment" :title="$t('experiment.index.config.title')">
    <div class="pl-6 w-full grid lg:grid-cols-2 lg:gap-10">
      <div class="pt-2">
        <p class="font-semibold pb-4">{{ $t('experiment.index.config.detail') }}</p>
        <SLTable :header="['Key', 'Value']" :data="getConfigs(experiment.config)" />
      </div>
      <div class="pt-2" v-if="summaries?.length !== 0">
        <p class="font-semibold pb-4">{{ $t('experiment.index.config.summarize') }}</p>
        <SLTable :header="['Key', 'Value']" :data="summaries" />
      </div>
      <div class="w-full min-h-30 flex justify-center items-center" v-else>
        <SLLoading />
      </div>
    </div>
    <div class="w-full h-6"></div>
  </ExtendBlock>
</template>

<script setup>
/**
 * @description: 实验的配置信息
 * @file: ExperimentConfig.vue
 * @since: 2023-12-11 17:07:31
 **/
import ExtendBlock from '@swanlab-vue/views/experiment/components/ExtendBlock.vue'
import SLTable from '@swanlab-vue/components/SLTable.vue'
import { inject } from 'vue'
import http from '@swanlab-vue/api/http'
import { ref } from 'vue'
import SLLoading from '@swanlab-vue/components/SLLoading.vue'

const experiment = inject('experiment')

// ---------------------------------- 转化实验配置为表格数据 ----------------------------------

const getConfigs = (config) => {
  const configs = []
  for (const key in config) {
    configs.push([key, config[key]])
  }
  return configs
}

// ---------------------------------- 获取实验的总结数据 ----------------------------------

const summaries = ref([])
const experimentId = inject('experimentId')
http.get(`/experiment/${experimentId.value}/summary`).then((res) => {
  summaries.value = res.data.summaries
  console.log(summaries.value)
})
</script>

<style lang="scss" scoped></style>
