import { defineStore } from 'pinia'
import { computed, ref } from 'vue'

export const useExperimentStroe = defineStore('charts', () => {
  /** state */
  // 当前实验
  const experiment = ref()
  /** getter */
  // 当前实验id
  const id = computed(() => experiment.value?.experiment_id)
  // 当前实验状态
  const status = computed(() => experiment.value?.status)
  // 当前实验颜色
  const color = computed(() => experiment.value?.color)
  // 默认颜色
  const defaultColor = computed(() => experiment.value?.default_color)

  /** action */
  // 设置当前实验状态
  const setStatus = (status) => {
    experiment.value.status = status
  }

  return {
    // state
    experiment,
    defaultColor,
    // getter
    id,
    status,
    color,
    // action
    setStatus
  }
})
