<template>
  <div class="w-full px-6 pt-6 text-dimmer">
    <!-- 实验标题 -->
    <div class="flex items-center">
      <span class="text-2xl font-semibold text-default pr-4">{{ experiment.name }}</span>
      <!-- <SLCopy :text="experiment.name" icon-class="w-5 h-5 text-dimmest cursor-pointer hover:text-dimmer mr-3" /> -->
      <!-- <SLIcon icon="setting" class="icon" /> -->
    </div>
    <!-- 实验描述 -->
    <div class="flex items-center pt-5" v-if="experiment?.description">
      <span>{{ experiment.description }}</span>
      <!-- <SLCopy
        :text="experiment.description"
        icon-class="w-5 h-5 text-dimmest cursor-pointer hover:text-dimmer ml-4 mr-3"
      /> -->
      <!-- <SLIcon icon="setting" class="icon" /> -->
    </div>
    <!-- 实验信息 -->
    <div class="flex justify-between pt-6 pb-2 flex-wrap">
      <!-- 实验相关 -->
      <div class="min-w-[400px] max-w-[50%] pr-4">
        <!-- 实验状态 -->
        <div class="flex pb-4">
          <div class="min-w-[150px]">{{ $t(`experiment.index.header.experiment_infos.status`) }}</div>
          <SLStatusLabel :name="experiment.name" :id="experiment.id" :status="experiment.status" />
          <!-- 停止按钮 -->
          <!-- <StopButton /> -->df3187f (feat: stop button ui)
        </div>
        <div v-for="item in experiment_infos" :key="item.title">
          <div class="flex pb-4" v-if="item.value">
            <div class="min-w-[150px]">{{ $t(`experiment.index.header.experiment_infos.${item.title}`) }}</div>
            <div v-if="!item.isLink" :title="item.value">{{ item.value }}</div>
            <a :href="item.value" target="_blank" class="hover:underline break-all" v-else>{{ item.value }}</a>
          </div>
        </div>
      </div>
      <!-- 系统相关 -->
      <div class="w-1/2 min-w-[400px]" v-if="experiment.system">
        <div v-for="item in experiment_device" :key="item.title" class="flex pb-4">
          <div class="min-w-[150px]">{{ $t(`experiment.index.header.experiment_device.${item.title}`) }}</div>
          <div :title="item.value">{{ item.value || 'Unkown' }}</div>
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
import SLCopy from '@swanlab-vue/components/SLCopy.vue'
import SLStatusLabel from '@swanlab-vue/components/SLStatusLabel.vue'
import { computed } from 'vue'
import { formatTime } from '@swanlab-vue/utils/time'
import { t } from '@swanlab-vue/i18n'
import { useExperimentStroe } from '@swanlab-vue/store'
import { ref } from 'vue'
import StopButton from './StopButton.vue'

const experiment = ref(useExperimentStroe().experiment)

// ---------------------------------- 实验信息 ----------------------------------

const experiment_infos = computed(() => {
  return [
    {
      title: 'start_time',
      value: formatTime(experiment.value.create_time)
    },
    {
      title: 'last_time',
      value: duration.value
    },
    {
      title: 'version',
      value: `v${experiment.value.version}`
    },
    {
      title: 'git',
      value: experiment.value.system.git_remote,
      isLink: true
    }
  ]
})

// ---------------------------------- 设备信息 ----------------------------------

const hardware = computed(() => {
  const prePath = 'experiment.index.header.experiment_device'
  const list = [
    experiment.value.system.cpu ? t(`${prePath}.cpu`, { value: experiment.value.system.cpu }) : '',
    experiment.value.system.gpu?.cores ? t(`${prePath}.gpu`, { value: experiment.value.system.gpu.cores }) : '',
    experiment.value.system.gpu?.type[0] ? t(`${prePath}.type`, { value: experiment.value.system.gpu.type[0] }) : ''
  ]
  return list.filter((item) => item !== '').join(' | ')
})
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
    },
    {
      title: 'executable',
      value: experiment.value.system.executable || ''
    },
    {
      title: 'hardware',
      value: hardware.value
    }
  ]
})

/**
 * 计算实验的持续时间
 */
const duration = computed(() => {
  const time1 = new Date(experiment.value.create_time)
  const currentTime = new Date()
  const time2 =
    experiment.value.status === 0
      ? new Date(currentTime.getTime() - 8 * 60 * 60 * 1000)
      : new Date(experiment.value.update_time)

  if (isNaN(time1.getTime()) || isNaN(time2.getTime())) {
    // 处理无效日期的情况
    return 'Invalid date'
  }

  const timeDifference = Math.abs(time2 - time1)

  const seconds = Math.floor(timeDifference / 1000) % 60
  const minutes = Math.floor(timeDifference / (1000 * 60)) % 60
  const hours = Math.floor(timeDifference / (1000 * 60 * 60)) % 24
  const days = Math.floor(timeDifference / (1000 * 60 * 60 * 24))

  const formattedTime = []

  if (days > 0) {
    formattedTime.push(`${days}${t('experiment.index.header.experiment_infos.time.day')}`)
  }

  if (hours > 0) {
    formattedTime.push(`${hours}${t('experiment.index.header.experiment_infos.time.hour')}`)
  }

  if (minutes > 0) {
    formattedTime.push(`${minutes}${t('experiment.index.header.experiment_infos.time.minute')}`)
  }

  if (seconds > 0) {
    formattedTime.push(`${seconds}${t('experiment.index.header.experiment_infos.time.second')}`)
  }

  return formattedTime.join('')
})
</script>

<style lang="scss" scoped>
.icon {
  @apply w-5 h-5 text-dimmest cursor-pointer hover:text-dimmer;
}
</style>
