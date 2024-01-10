<template>
  <div class="w-full flex gap-5 justify-between">
    <SLSearch @input="input" class="max-w-[400px]" />
    <div class="flex gap-3">
      <SLButton hollow @click="copy">{{ $t('experiment.func-bar.copy') }}</SLButton>
      <SLButton hollow @click="download">{{ $t('experiment.func-bar.download') }}</SLButton>
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
import { copyTextToClipboard, downloadStringAsFile } from '@swanlab-vue/utils/browser'

const props = defineProps({
  // 复制、下载的内容
  content: {
    type: String,
    default: ''
  },
  // 下载时使用的文件名
  filename: {
    type: String,
    default: 'swanlab.file'
  }
})

const emits = defineEmits(['input', 'download'])

// 输入触发
const input = (value) => {
  emits('input', value)
}

// 复制
const copy = () => {
  copyTextToClipboard(props.content)
  message.success(t('experiment.func-bar.copy-success'))
}

// 下载
const download = () => {
  downloadStringAsFile(props.content, props.filename)
}
</script>

<style lang="scss" scoped>
button {
  @apply rounded-lg px-3 py-1;
}
</style>
