<template>
  <ExtendBlock class="pt-2 pr-10" icon="config" :title="$t('experiment.index.config.title')">
    <div class="pl-6 w-full grid lg:grid-cols-2 lg:gap-10">
      <div class="pt-2">
        <div class="flex items-center pb-4">
          <p class="font-semibold pr-2">{{ $t('experiment.index.config.detail') }}</p>
          <SLHelp>这是初始化实验时的初始化配置</SLHelp>
        </div>
        <SLTable class="max-w-[600px]" :header="['Key', 'Value']" :data="getConfigs(experiment.config)" />
      </div>
      <div class="pt-2" v-if="summaries?.length !== 0">
        <div class="flex items-center pb-4">
          <p class="font-semibold pr-2">{{ $t('experiment.index.config.summarize') }}</p>
          <SLHelp>这是每个tag中最后一个step的数据</SLHelp>
        </div>
        <SLTable class="max-w-[600px]" :header="['Key', 'Value']" :data="summaries" />
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
import SLHelp from '@swanlab-vue/components/SLHelp.vue'
import http from '@swanlab-vue/api/http'
import { ref } from 'vue'
import SLLoading from '@swanlab-vue/components/SLLoading.vue'
import { useExperimentStroe } from '@swanlab-vue/store'
import { inject } from 'vue'

const experiment = ref(useExperimentStroe().experiment)

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
const experimentId = ref(useExperimentStroe().id)
const show_error = inject('show_error')
http
  .get(`/experiment/${experimentId.value}/summary`)
  .then((res) => {
    summaries.value = res.data.summaries
    console.log(summaries.value)
  })
  .catch((error) => {
    console.error(error)
    show_error(error.data?.code || 500)
  })
</script>

<style lang="scss" scoped></style>
