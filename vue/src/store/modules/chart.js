import { defineStore } from 'pinia'
import { ref } from 'vue'

export const useChartsStore = defineStore('charts', () => {
  /** state */
  const charts = ref([])
  /** getter */

  /** action */

  return {
    charts
  }
})
