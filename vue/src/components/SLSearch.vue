<template>
  <div class="search" :class="focused ? ' border-primary-default' : 'hover:border-primary-dimmer'">
    <SLIcon class="w-4 h-4 shrink-0" icon="search"></SLIcon>
    <input
      type="text"
      v-model="value"
      class="w-full outline-none text-sm text-dimmer"
      :placeholder="placeholder"
      @input="input"
      @keydown.enter="$emit('search', value)"
      @focus="focused = true"
      @blur="focused = false"
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
  }
})

const value = ref('')
const focused = ref(false)

const input = debounce(() => {
  emits('input', value.value)
  emits('update:modelValue', value.value)
}, props.dealy)
</script>

<style lang="scss" scoped>
.search {
  @apply w-full border flex gap-2 justify-between items-center py-2 px-2 rounded-lg transition-all;
}
</style>
