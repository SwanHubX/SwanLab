import { defineStore } from 'pinia'
import { computed, ref } from 'vue'

export const useExperimentStroe = defineStore('charts', () => {
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

  /** action */
  // 设置当前实验状态
  const setStatus = (status) => {
    experiment.value.status = status
  }
  // 修改更新时间
  const setUpateTIme = (time) => {
    experiment.value.update_time = time
  }

  // 修改实验信息
  const setExperiment = (x) => {
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
    // action
    setStatus,
    setUpateTIme,
    setExperiment
  }
})
