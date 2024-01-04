import { defineStore } from 'pinia'
import { ref, computed } from 'vue'

export const useProjectStore = defineStore('project', () => {
  /** state */
  const project = ref()
  /** getter */
  const name = computed(() => (project.value?.name ? project.value.name : 'my-machine-learning-project'))
  const description = computed(() => (project.value?.description ? project.value.description : ''))
  const experiments = computed(() => project.value?.experiments)
  const sum = computed(() => project.value?.experiments.length)
  const createTime = computed(() => project.value?.create_time)
  const updateTime = computed(() => project.value?.update_time)

  /** action */
  const setProject = (p) => {
    p.experiments.reverse()
    project.value = p
  }
  /**
   * 清空 project
   */
  const clearProject = () => {
    project.value = null
  }
  /**
   * 设置指定实验的状态
   * @param { number } id 实验id
   * @param { string } status 实验状态
   */
  const setExperimentStatus = (id, status) => {
    // 不需要绝对等于
    const experiment = experiments.value.find((e) => e.experiment_id == id)
    experiment.status = status
  }
  /**
   * 重置实验列表中某一个实验的信息
   * @param {number} id 燕燕唯一id
   * @param {object} newInfo 涵盖字段为experiment中字段的子集，不需要修改的不必传
   * @param {boolean} overwirte 设置true后，newInfo直接覆盖experiment，请保证newInfo的完整性
   * @returns
   */
  const setExperimentInfo = (id, newInfo, overwirte) => {
    let experiment = experiments.value.find((e) => e.experiment_id == id)
    if (overwirte) {
      experiment = newInfo
      return
    }
    // 遍历 newInfo 中的 key
    for (const key in newInfo) {
      // 如果experiment中没有这个key，说明key不合法
      if (!experiment[key]) continue
      // key 合法，更新该字段
      experiment[key] = newInfo[key]
    }
  }

  return {
    sum,
    name,
    description,
    experiments,
    createTime,
    updateTime,
    setProject,
    clearProject,
    setExperimentStatus,
    setExperimentInfo
  }
})
