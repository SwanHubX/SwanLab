<template>
  <svg aria-hidden="true">
    <use :xlink:href="'#icon-' + icon"></use>
  </svg>
</template>

<script setup>
/**
 * @description: 对iconfont图标的使用进行封装
 * @file: SLIcon.vue
 * @since: 2023-12-04 16:07:51
 **/
import '@swanlab-vue/assets/iconfont.js'
defineProps({
  icon: {
    type: String,
    required: true
  }
})
</script>
