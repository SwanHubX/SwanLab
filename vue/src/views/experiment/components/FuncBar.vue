<template>
  <div class="w-full flex gap-5 justify-between">
    <SLSearch @input="input" class="max-w-[400px]" />
    <div class="flex gap-3">
      <SLButton hollow @click="copy">{{ $t('experiment.func-bar.copy') }}</SLButton>
      <SLButton hollow>{{ $t('experiment.func-bar.download') }}</SLButton>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 功能条：搜索、复制、下载
 * @file: FuncBar.vue
 * @since: 2024-01-10 12:54:19
 **/

import SLSearch from '@swanlab-vue/components/SLSearch.vue'
import { message } from '@swanlab-vue/components/message'
import { t } from '@swanlab-vue/i18n'
import { copyTextToClipboard } from '@swanlab-vue/utils/browser'

const props = defineProps({
  copyText: {
    type: String,
    default: ''
  },
  downloadName: {
    type: String,
    default: 'swanlab.file'
  },
  downloadContent: {
    type: String,
    default: 'SwanLab'
  }
})

const emits = defineEmits(['input', 'download'])

const input = (value) => {
  emits('input', value)
}

const copy = () => {
  copyTextToClipboard(props.copyText)
  message.success(t('experiment.func-bar.copy-success'))
}

const download = () => {}
</script>

<style lang="scss" scoped>
button {
  @apply rounded-lg px-3 py-1;
}
</style>
