<template>
  <a class="sl-status-label" :class="config.class" :href="url" @click.prevent="handleClick">
    {{ config.label }}
  </a>
</template>

<script setup>
/**
 * @description: 实验状态展示，status为0代表正在运行，为1代表运行成功，-1代表运行失败
 * @param { Number } status 实验状态
 * @param { Number } id 实验id，如果不传入id，点击状态标签不会跳转
 * @file: StatusLabel.vue
 * @since: 2023-12-08 20:24:49
 **/

import { t } from '@swanlab-vue/i18n'
import { computed } from 'vue'
import { useRouter } from 'vue-router'
const router = useRouter()
const props = defineProps({
  status: {
    type: Number,
    required: true
  },
  id: {
    default: null
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

// ---------------------------------- 将id转译为url ----------------------------------

const url = computed(() => {
  if (!props.id) return undefined
  return '/experiment/' + props.id
})

// ---------------------------------- 代理a标签点击事件 ----------------------------------
const handleClick = () => {
  if (!url.value) return
  router.push(url.value)
}
</script>

<style lang="scss" scoped>
.sl-status-label {
  @apply px-3 py-1 rounded-full text-sm text-center;
}
.stoped {
  @apply bg-negative-dimmest text-negative-default;
}

.finished {
  @apply bg-positive-dimmest text-positive-higher;
}

.running {
  @apply bg-primary-dimmest text-primary-higher;
}
</style>
