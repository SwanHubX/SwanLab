<template>
  <div class="w-full flex justify-between gap-5">
    <div class="flex items-center w-full">
      <SLSearch @input="input" :placeholder="placeholder" class="max-w-[400px]" />
      <div class="flex items-center px-4" v-show="targets?.length">
        <span>{{ _modelValue }} / {{ targets?.length }}</span>
        <div class="w-3 flex-shrink-0 flex-col flex ml-2">
          <SLIcon icon="down" class="w-full h-3 -rotate-180 -mb-1" @click="_modelValue++" />
          <SLIcon icon="down" class="w-full aspect-square" @click="_modelValue--" />
        </div>
      </div>
    </div>
    <div class="flex gap-3">
      <SLButton hollow @click="copy">
        <SLIcon icon="copy" class="icon"></SLIcon>
        {{ $t('experiment.func-bar.copy') }}
      </SLButton>
      <SLButton hollow @click="download">
        <SLIcon icon="download" class="icon"></SLIcon>
        {{ $t('experiment.func-bar.download') }}
      </SLButton>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 功能条：搜索、复制、下载
 * @file: FuncBar.vue
 * @since: 2024-01-10 12:54:19
 **/

import { computed } from 'vue'
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
  },
  // 搜索栏的占位
  placeholder: {
    type: String,
    default: 'Search...'
  },
  // 目标行
  targets: {
    type: Array,
    default: () => []
  },
  // 当前选中的搜索目标
  modelValue: {
    type: Number,
    default: 1
  }
})

const emits = defineEmits(['input', 'download', 'update:modelValue'])

const _modelValue = computed({
  get() {
    return props.modelValue
  },
  set(v) {
    if (v < 1) {
      emits('update:modelValue', 1)
    } else if (v > props.targets.length) {
      emits('update:modelValue', props.targets.length)
    } else {
      emits('update:modelValue', v)
    }
  }
})

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
  @apply rounded-lg px-3 py-1 flex items-center gap-2;
}

.icon {
  @apply w-4 h-4;
}
</style>
