<template>
  <MenuItem :as="as" class="p-3 hover:bg-higher cursor-pointer">
    <slot></slot>
  </MenuItem>
</template>

<script setup>
/**
 * SLMenuItem - 菜单项组件，封装headlessui的MenuItem组件
 */
import { MenuItem } from '@headlessui/vue'

defineProps({
  as: {
    type: String,
    default: 'div'
  }
})
</script>
