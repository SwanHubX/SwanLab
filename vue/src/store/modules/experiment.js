import { defineStore } from 'pinia'
import { computed, ref } from 'vue'
import { t } from '@swanlab-vue/i18n'

export const useExperimentStroe = defineStore('experiment', () => {
  /** state */
  // 当前实验
  const experiment = ref()
  // 当前实验的charts
  const charts = ref()
  /** getter */
  // 当前实验id
  const id = computed(() => experiment.value?.experiment_id)
  // 实验名称
  const name = computed(() => experiment.value?.name)
  // 实验描述
  const description = computed(() => experiment.value?.description)
  // 当前实验状态
  const status = computed(() => experiment.value?.status)
  // 当前实验颜色
  const color = computed(() => experiment.value?.color)
  // 默认颜色
  const defaultColor = computed(() => experiment.value?.default_color)
  // 是否running
  const isRunning = computed(() => status.value === 0)
  // 持续时间
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

    return formattedTime.join('') || 'less than 1s'
  })

  /** action */
  // 设置当前实验状态
  const setStatus = (status) => {
    experiment.value.status = status
  }
  // 修改更新时间
  const setUpateTime = (time) => {
    experiment.value.update_time = time
  }

  // 修改实验信息
  const setExperiment = (x) => {
    x.default_color = defaultColor.value
    experiment.value = x
  }

  return {
    // state
    experiment,
    charts,
    // getter
    id,
    name,
    description,
    status,
    color,
    defaultColor,
    isRunning,
    duration,
    // action
    setStatus,
    setUpateTime,
    setExperiment
  }
})
