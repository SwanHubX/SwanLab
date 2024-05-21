<template>
  <div class="flex items-center gap-2 cursor-pointer" @click="!disabled && $emit('update:checked', !checked)">
    <input type="checkbox" class="border cursor-pointer w-4 h-4" :checked="checked && !disabled" :disabled="disabled" />
    <span class="select-none truncate">{{ label }}</span>
  </div>
</template>

<script setup>
/**
 * @description: 选项框
 * @file: SLCheck.vue
 * @since: 2024-01-16 18:13:34
 **/

defineProps({
  label: {
    type: String,
    default: ''
  },
  checked: {
    type: Boolean,
    default: false
  },
  disabled: {
    type: Boolean,
    default: false
  }
})

defineEmits(['update:checked'])
</script>

<style lang="scss" scoped>
input[type='checkbox']:disabled {
  @apply cursor-not-allowed;
}
</style>
