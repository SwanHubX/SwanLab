<template>
  <div class="w-full border flex gap-2 justify-between items-center py-1 px-2 rounded-lg">
    <SLIcon class="w-4 h-4" icon="search"></SLIcon>
    <input
      type="text"
      v-model="value"
      class="w-full outline-none"
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

const emits = defineEmits(['input', 'search'])

const props = defineProps({
  dealy: {
    type: String,
    default: '500'
  }
})

const value = ref('')

const input = debounce(() => {
  emits('input', value.value)
}, props.dealy)
</script>

<style lang="scss" scoped></style>
