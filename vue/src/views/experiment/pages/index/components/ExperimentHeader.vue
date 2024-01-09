<template>
  <div class="w-full px-6 pb-2 text-default relative">
    <!-- 实验信息 -->
    <div class="flex justify-between pt-6 pb-2 flex-wrap">
      <!-- 实验相关 -->
      <div class="min-w-[400px] max-w-[50%] pr-4">
        <!-- 实验状态 -->
        <div class="flex items-center pb-4">
          <div class="title">{{ $t(`experiment.index.header.experiment_infos.status`) }}</div>
          <SLStatusLabel :name="experiment.name" :id="experiment.id" :status="experiment.status" />
          <!-- 停止按钮 -->
          <StopButton v-if="experimentStore.isRunning" />
        </div>
        <div v-for="item in experiment_infos" :key="item.title">
          <div class="flex pb-4" v-if="item.value || item.title === 'git'">
            <div class="title">{{ $t(`experiment.index.header.experiment_infos.${item.title}`) }}</div>
            <div v-if="!item.isLink" :title="item.value">{{ item.value }}</div>
            <div v-else-if="item.title === 'git' && !item.value" class="text-dimmest">None</div>
            <a :href="item.value" target="_blank" class="hover:underline break-all" v-else>{{ item.value }}</a>
          </div>
        </div>
      </div>
      <!-- 系统相关 -->
      <div class="w-1/2 min-w-[400px]" v-if="experiment.system">
        <div v-for="item in experiment_device" :key="item.title" class="flex pb-4">
          <div class="title">{{ $t(`experiment.index.header.experiment_device.${item.title}`) }}</div>
          <div :title="item.value">{{ item.value || 'Unkown' }}</div>
        </div>
        <SLButton theme="primary" hollow class="rounded-lg">
          <router-link to="env" class="block px-3 py-1">All Details</router-link>
        </SLButton>
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
import SLStatusLabel from '@swanlab-vue/components/SLStatusLabel.vue'
import { computed, ref } from 'vue'
import { formatTime } from '@swanlab-vue/utils/time'
import { t } from '@swanlab-vue/i18n'
import { useExperimentStroe, useProjectStore } from '@swanlab-vue/store'
import http from '@swanlab-vue/api/http'
import { useRouter } from 'vue-router'
import { inject } from 'vue'
import { message } from '@swanlab-vue/components/message'
import StopButton from './StopButton.vue'

const experimentStore = useExperimentStroe()
const experiment = ref(experimentStore.experiment)

const projectStore = useProjectStore()

const showErrorView = inject('showErrorView')

const router = useRouter()

// ---------------------------------- 实验信息 ----------------------------------

const experiment_infos = computed(() => {
  return [
    {
      title: 'start_time',
      value: formatTime(experiment.value.create_time)
    },
    {
      title: 'last_time',
      value: experimentStore.duration
    },
    {
      title: 'version',
      value: `v${experiment.value.version}`
    }
  ]
})

// ---------------------------------- 设备信息 ----------------------------------

// const hardware = computed(() => {
//   const prePath = 'experiment.index.header.experiment_device'
//   const list = [
//     experiment.value.system.cpu ? t(`${prePath}.cpu`, { value: experiment.value.system.cpu }) : '',
//     experiment.value.system.gpu?.cores ? t(`${prePath}.gpu`, { value: experiment.value.system.gpu.cores }) : '',
//     experiment.value.system.gpu?.type[0] ? t(`${prePath}.type`, { value: experiment.value.system.gpu.type[0] }) : ''
//   ]
//   return list.filter((item) => item !== '').join(' | ')
// })

const experiment_device = computed(() => {
  return [
    {
      title: 'hostname',
      value: experiment.value.system.hostname || ''
    },
    {
      title: 'os',
      value: experiment.value.system.os || ''
    },
    {
      title: 'python',
      value: experiment.value.system.python || ''
    }
  ]
})
</script>

<style lang="scss" scoped>
.icon {
  @apply w-5 h-5 text-dimmest cursor-pointer hover:text-dimmer;
}

.title {
  @apply min-w-[150px] text-dimmest;
}
</style>
