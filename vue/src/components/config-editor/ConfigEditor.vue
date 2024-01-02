<template>
  <button
    class="flex gap-2 items-center px-2 py-1 border rounded-lg transition-all hover:bg-primary-default hover:text-white"
    @click="() => (showModal = true)"
  >
    <SLIcon icon="modify" class="w-4 h-5" />
    {{ $t('common.config-editor.button') }}
  </button>
  <SLModal class="px-10 py-6" max-w="550" v-model="showModal">
    <!-- 如果后续项目和实验的可修改内容发生改变，这样容易重构一些 -->
    <EditorWrap :type="type" @confirm="confirm"></EditorWrap>
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
import SLIcon from '../SLIcon.vue'

const showModal = ref(false)

defineProps({
  type: {
    type: String,
    validator: (value) => {
      return ['project', 'experiment'].includes(value)
    }
  }
})

const emits = defineEmits(['modify'])

const confirm = async (newV) => {
  new Promise((resolve) => {
    emits('modify', newV, () => {
      showModal.value = false
      resolve()
    })
  })
}
</script>

<style lang="scss" scoped></style>
