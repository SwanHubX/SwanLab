import { defineStore } from 'pinia'
import { ref, computed } from 'vue'

export const useProjectStore = defineStore('project', () => {
  /** state */
  const project = ref()
  /** getter */
  const name = computed(() => 'my-machine-learning-project')
  const description = computed(() => '实验描述也没有配置')
  const experiments = computed(() => project.value?.experiments)
  const sum = computed(() => project.value?._sum)
  const createTime = computed(() => project.value?.create_time)
  const updateTime = computed(() => project.value?.update_time)

  /** action */
  const setProject = (p) => {
    p.experiments.reverse()
    project.value = p
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

  return {
    sum,
    name,
    description,
    experiments,
    createTime,
    updateTime,
    setProject,
    setExperimentStatus
  }
})
