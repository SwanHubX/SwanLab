<template>
  <div class="swan-search" :class="{ 'icon-reverse': reverse }">
    <SLIcon class="swan-icon" icon="search" />
    <input
      type="text"
      v-model="value"
      class="focus:border-primary-default hover:border-primary-dimmer"
      :placeholder="placeholder"
      @input="input"
      @keydown.enter="$emit('search', value)"
    />
  </div>
</template>

<script setup>
/**
 * @description: 搜索
 * @file: SLSearch.vue
 * @since: 2024-01-08 16:36:14
 **/

import { ref } from 'vue'
import SLIcon from './SLIcon.vue'
import { debounce } from '@swanlab-vue/utils/common'
import { t } from '@swanlab-vue/i18n'

const emits = defineEmits(['input', 'search', 'update:modelValue'])

const props = defineProps({
  dealy: {
    type: String,
    default: '0'
  },
  placeholder: {
    type: String,
    default: t('common.search.placeholder')
  },
  modelValue: {
    type: String,
    default: ''
  },
  reverse: {
    type: Boolean,
    default: false
  }
})

const value = ref('')

const input = debounce(() => {
  emits('input', value.value)
  emits('update:modelValue', value.value)
}, props.dealy)
</script>

<style lang="scss" scoped>
.swan-search {
  @apply relative w-full;
  input {
    @apply w-full outline-none text-sm text-dimmer border p-2 rounded-lg transition-all;
    @apply pl-7;
  }

  .swan-icon {
    @apply w-4 h-4 shrink-0;
    @apply absolute top-1/2 left-2 transform -translate-y-1/2;
  }
}

.icon-reverse {
  input {
    @apply pl-2 pr-7;
  }
  .swan-icon {
    @apply right-2 left-auto;
  }
}
</style>
