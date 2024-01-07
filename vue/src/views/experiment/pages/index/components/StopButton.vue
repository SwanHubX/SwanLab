<template>
  <button
    class="p-1.5 rounded-full border-2 transition-all duration-200 ml-2"
    :class="showColor ? 'border-negative-default' : 'cursor-not-allowed'"
    @mouseover="() => (hover = true)"
    @mouseout="() => (hover = false)"
    @click="stop"
  >
    <div
      class="w-2.5 h-2.5 transition-all duration-200"
      :class="showColor ? 'bg-negative-default' : 'bg-overlay opacity-50'"
    ></div>
  </button>
</template>

<script setup>
/**
 * @description: 停止实验按钮
 * @file: StopButton.vue
 * @since: 2023-12-30 20:30:30
 **/
import { useExperimentStroe, useProjectStore } from '@swanlab-vue/store'
import { ref } from 'vue'
import { computed } from 'vue'
import { confirm } from '@swanlab-vue/components/comfirm'
import http from '@swanlab-vue/api/http'
import { t } from '@swanlab-vue/i18n'

// ---------------------------------- 弹窗相关 ----------------------------------

const experiment = useExperimentStroe()
const id = experiment.id
const status = experiment.status

const hover = ref(false) // 是否hover
// 展示hover时可点击样式
const showColor = computed(() => {
  return status === 0 && hover.value
})

// ---------------------------------- 确认删除 ----------------------------------

const stop = () => {
  confirm(
    t('experiment.index.header.stop.modal.title'),
    t('experiment.index.header.stop.modal.text'),
    t('experiment.index.header.stop.button')
  ).then(stop_experiment)
}

const stop_experiment = async () => {
  const { data } = await http.get(`/experiment/${id}/stop`)
  if (!data) return
  experiment.setUpateTIme(data.update_time)
  useProjectStore().setExperimentStatus(id, -1)
}
</script>

<style lang="scss" scoped></style>
