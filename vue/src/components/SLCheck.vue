<template>
  <div class="flex items-center gap-2 cursor-pointer" @click="() => (checked = !checked)">
    <input type="checkbox" class="border w-4 h-4" :checked="checked" />
    <span class="select-none truncate">{{ label }}</span>
  </div>
</template>

<script setup>
/**
 * @description: 选项框
 * @file: SLCheck.vue
 * @since: 2024-01-16 18:13:34
 **/
import { ref } from 'vue'

defineProps({
  label: {
    type: String,
    default: ''
  }
})

const checked = ref(false)
</script>

<style lang="scss" scoped></style>
