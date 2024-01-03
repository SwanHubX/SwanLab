<template>
  <SLModal v-model="show"> </SLModal>
</template>

<script setup>
import { ref } from 'vue'
import SLModal from '../SLModal.vue'

/**
 * @description: 检查组件，用于用户进行危险操作时的确认，基于SLModal.vue实现弹窗效果
 * @file: SLCheck.vue
 * @since: 2024-01-03 14:42:37
 **/
defineProps({
  title: {
    type: String,
    default: '提示'
  },
  text: {
    type: String,
    default: '确认删除吗？'
  },
  config: {
    buttonText: {
      type: Object,
      default: {}
    }
  }
})

const show = ref(true)
</script>

<style lang="scss" scoped></style>
