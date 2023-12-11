<template>
  <div class="w-full px-6 pt-6 text-dimmer">
    <!-- 实验标题 -->
    <div class="flex items-center">
      <span class="text-2xl font-semibold text-default">{{ experiment.name }}</span>
      <SLStatusLabel :name="experiment.name" :id="experiment.id" :status="experiment.status" class="mx-4" />
      <SLIcon icon="copy" class="icon mr-3" @click="copyTextToClipboard(experiment.name)" />
      <SLIcon icon="setting" class="icon" />
    </div>
    <!-- 实验描述 -->
    <div class="flex items-center pt-5">
      <span>{{ experiment.description }}</span>
      <SLIcon icon="copy" class="icon ml-4 mr-3" @click="copyTextToClipboard(experiment.description)" />
      <SLIcon icon="setting" class="icon" />
    </div>
    <!-- 实验信息 -->
    <div class="flex justify-between pt-6 pb-2 flex-wrap">
      <!-- 实验相关 -->
      <div class="w-1/2 min-w-[400px]">
        <div v-for="item in experiment_infos" :key="item.title" class="flex pb-4">
          <div class="min-w-[150px]">{{ $t(`experiment.header.experiment_infos.${item.title}`) }}</div>
          <div class="">{{ item.value }}</div>
        </div>
      </div>
      <!-- 系统相关 -->
      <div class="w-1/2 min-w-[400px]">
        <div v-for="item in experiment_device" :key="item.title" class="flex pb-4">
          <div class="min-w-[180px]">{{ $t(`experiment.header.experiment_device.${item.title}`) }}</div>
          <div class="">{{ item.value }}</div>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 实验概览页-信息头部
 * @file: ExperimentHeader.vue
 * @since: 2023-12-11 14:43:51
 **/
import { copyTextToClipboard } from '@swanlab-vue/utils/browser'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import SLStatusLabel from '@swanlab-vue/components/SLStatusLabel.vue'
import { computed } from 'vue'
import { convertUtcToLocal, transTime } from '@swanlab-vue/utils/time'
import { inject } from 'vue'

const experiment = inject('experiment')
console.log(experiment.value)

const experiment_infos = computed(() => {
  return [
    {
      title: 'start_time',
      value: transTime(convertUtcToLocal(experiment.value.create_time))
    },
    {
      title: 'last_time',
      value: experiment.value.update_time - experiment.value.create_time
    }
  ]
})

const experiment_device = computed(() => {
  return [
    {
      title: 'hostname',
      value: experiment.value.hostname || '未知'
    },
    {
      title: 'os',
      value: experiment.value.os || '未知'
    },
    {
      title: 'python',
      value: experiment.value.python || '未知'
    }
  ]
})
</script>

<style lang="scss" scoped>
.icon {
  @apply w-5 h-5 text-dimmest cursor-pointer hover:text-dimmer;
}
</style>
