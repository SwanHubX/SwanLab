<template>
  <div class="flex items-center gap-2">
    <input type="checkbox" class="border" :value="value" />
    <span>{{ label }}</span>
  </div>
</template>

<script setup>
/**
 * @description: 选项框
 * @file: SLCheck.vue
 * @since: 2024-01-16 18:13:34
 **/
import { ref } from 'vue'

const value = ref(false)

defineProps({
  label: {
    type: String,
    default: ''
  }
})
</script>

<style lang="scss" scoped></style>
