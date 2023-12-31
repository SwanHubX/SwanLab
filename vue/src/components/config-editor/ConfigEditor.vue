<template>
  <button class="px-2 py-1 border rounded-lg" @click="() => (showModal = true)">Edit</button>
  <SLModal class="px-10 py-6" max-w="550" v-model="showModal">
    <!-- 如果后续项目和实验的可修改内容发生改变，这样容易重构一些 -->
    <EditorWrap :type="type"></EditorWrap>
  </SLModal>
</template>

<script setup>
/**
 * @description: 编辑项目/实验信息
 * @file: InfoEditor.vue
 * @since: 2023-12-30 23:55:02
 **/
import { ref } from 'vue'
import SLModal from '../SLModal.vue'
import EditorWrap from './EditorWrap.vue'

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
