import { defineStore } from 'pinia'
import { computed, ref } from 'vue'
import { getDuration } from '@swanlab-vue/utils/time'

export const useExperimentStore = defineStore('experiment', () => {
  /** state */
  // 当前实验
  const experiment = ref()
  // 当前实验的charts
  const charts = ref()
  /** getter */
  // 当前实验id
  const id = computed(() => experiment.value?.id)
  // 实验名称
  const name = computed(() => experiment.value?.name)
  // 实验描述
  const description = computed(() => experiment.value?.description)
  // 当前实验状态
  const status = computed(() => experiment.value?.status)
  // 当前实验颜色
  const color = computed(() => experiment.value?.light)
  // 是否running
  const isRunning = computed(() => status.value === 0)
  // 持续时间
  const duration = computed(() => {
    return getDuration(experiment.value) || ''
  })

  /** action */
  // 设置当前实验状态
  const setStatus = (status) => {
    experiment.value.status = status
  }
  // 修改完成时间
  const setFinishTime = (time) => {
    experiment.value.finish_time = time
  }
  // 修改实验信息 => 全部信息
  const setExperiment = (x) => {
    experiment.value = x
  }
  // 修改实验信息 => 部分信息
  const setExperimentPartial = (data) => {
    for (const key in data) {
      // 如果experiment中没有这个key，说明key不合法
      if (experiment.value[key] == undefined) continue
      // key 合法，更新该字段
      experiment.value[key] = data[key]
    }
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
    isRunning,
    duration,
    // action
    setStatus,
    setFinishTime,
    setExperiment,
    setExperimentPartial
  }
})
