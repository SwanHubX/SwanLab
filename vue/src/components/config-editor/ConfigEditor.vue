<template>
  <button class="px-2 py-1 border rounded-lg" @click="() => (showModal = true)">Edit</button>
  <SLModal class="px-10 pt-10 pb-5" max-w="550" v-model="showModal"></SLModal>
</template>

<script setup>
/**
 * @description: 编辑项目/实验信息
 * @file: InfoEditor.vue
 * @since: 2023-12-30 23:55:02
 **/
import { ref } from 'vue'
import SLModal from '../SLModal.vue'

const showModal = ref(false)

defineProps({
  type: {
    type: String,
    validator: (value) => {
      return ['project', 'experiment'].includes(value)
    }
  }
})
</script>

<style lang="scss" scoped></style>
