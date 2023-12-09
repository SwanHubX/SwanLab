<template>
  <RouterLink class="px-3 py-1 rounded-full text-sm" :class="config.class" :to="`/experiment/${id}`">
    {{ config.label }}
  </RouterLink>
</template>

<script setup>
/**
 * @description: 实验状态展示，status为0代表正在运行，为1代表运行成功，-1代表运行失败
 * @file: StatusLabel.vue
 * @since: 2023-12-08 20:24:49
 **/

import { t } from '@swanlab-vue/i18n'
import { computed } from 'vue'
const props = defineProps({
  name: {
    type: String,
    required: true
  },
  status: {
    type: Number,
    required: true
  },
  id: {
    required: true
  }
})

// ---------------------------------- 将status转译为对应的状态提示和class ----------------------------------

const config = computed(() => {
  const status = t('experiment.status.' + props.status)
  let className = ''
  switch (props.status) {
    case 0:
      className = 'running'
      break
    case 1:
      className = 'finished'
      break
    case -1:
      className = 'stoped'
      break
  }
  if (!className) throw new Error('Invalid status:' + props.status + 'in StatusLabel.vue')

  return {
    label: status,
    class: className
  }
})
</script>

<style lang="scss" scoped>
.stoped {
  @apply bg-negative-dimmest text-negative-default;
}

.finished {
  @apply bg-positive-higher text-positive-dimmer;
}

.running {
  @apply bg-primary-dimmest text-primary-default;
}
</style>
