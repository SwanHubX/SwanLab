<template>
  <button class="underline text-dimmest ml-2" :class="showColor ? '' : 'cursor-not-allowed'" @click="stop">
    is over?
  </button>
</template>

<script setup>
/**
 * @description: 停止实验按钮
 * @file: StopButton.vue
 * @since: 2023-12-30 20:30:30
 **/
import { useExperimentStore, useProjectStore } from '@swanlab-vue/store'
import { computed } from 'vue'
import { confirm } from '@swanlab-vue/components/confirm'
import http from '@swanlab-vue/api/http'
import { t } from '@swanlab-vue/i18n'
import { message } from '@swanlab-vue/components/message'

// ---------------------------------- 弹窗相关 ----------------------------------

const experiment = useExperimentStore()
const id = experiment.id
const status = experiment.status

// 展示hover时可点击样式
const showColor = computed(() => {
  return status === 0
})

// ---------------------------------- 确认停止 ----------------------------------

const stop = () => {
  confirm(t('experiment.index.header.stop.modal.title'), t('experiment.index.header.stop.modal.text')).then(
    stop_experiment
  )
}

const stop_experiment = async () => {
  const { data } = await http.get(`/experiment/${id}/stop`)
  if (!data) return
  experiment.setStatus(data.status)
  experiment.setFinishTime(data.finish_time)
  useProjectStore().setExperimentStatus(id, -1, data.finish_time)
  message.success('State changed')
}
</script>

<style lang="scss" scoped></style>
