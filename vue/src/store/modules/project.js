import { defineStore } from 'pinia'
import { ref, computed } from 'vue'

export const useProjectStore = defineStore('project', () => {
  /** state */
  const project = ref()
  /** getter */
  const experiments = computed(() => project.value?.experiments)
  const sum = computed(() => project.value?._sum)
  const createTime = computed(() => project.value?.create_time)
  const updateTime = computed(() => project.value?.update_time)

  /** action */
  const setProject = (p) => {
    project.value = p
  }

  return {
    sum,
    experiments,
    createTime,
    updateTime,
    setProject
  }
})
