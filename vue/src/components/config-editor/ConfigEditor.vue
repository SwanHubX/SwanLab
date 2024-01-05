<template>
  <SLButton
    class="flex gap-2 items-center px-3 py-1.5 border rounded-lg text-sm"
    @click="() => (showModal = true)"
    theme="primary"
    hollow
    :class="disabled ? 'cursor-not-allowed' : ''"
    :disabled="disabled"
    :disabled-tip="$t('common.config-editor.not-allowed')"
  >
    <SLIcon icon="modify" class="w-3.5 h-3.5" />
    {{ $t('common.config-editor.button') }}
  </SLButton>
  <SLModal class="px-10 py-6" max-w="550" v-model="showModal">
    <!-- 如果后续项目和实验的可修改内容发生改变，这样容易重构一些 -->
    <EditorWrap :type="type" @confirm="confirm"></EditorWrap>
  </SLModal>
</template>

<script setup>
/**
 * @description: 编辑项目/实验信息
 * 目前文字部分只适配了项目和实验，如果之后需要使用这套弹窗模板修改别的部分的信息
 * 只需要增加i18n的文字适配，点击后的操作使用 modify 这个 emit 暴露了出去
 * @file: InfoEditor.vue
 * @since: 2023-12-30 23:55:02
 **/
import { ref } from 'vue'
import SLModal from '../SLModal.vue'
import EditorWrap from './EditorWrap.vue'
import SLIcon from '../SLIcon.vue'
import SLButton from '../SLButton.vue'

// 是否展示弹窗
const showModal = ref(false)

defineProps({
  type: {
    type: String,
    validator: (value) => {
      // 目前type只适配了project和experiment两种模式
      return ['project', 'experiment'].includes(value)
    }
  },
  disabled: {
    type: Boolean,
    default: false
  }
})

const emits = defineEmits(['modify'])

/**
 * 点击修改后触发的确认函数
 * @param {obj} newV 需要修改的值
 * @param.name {string}
 * @param.description {string}
 *
 * 这里使用了promise，是因为 emit 触发的 modify 很可能是异步的，而 emit 触发是同步的
 * 该组件在点击确认按钮后会有两个缓冲操作：
 * 1. 按钮显示处理中
 * 2. 弹窗并不立即消失
 * 可以看到，触发 modify 的时候会传递两个参数，一个是修改后的值，一个是回调函数
 * 只有在外部的处理函数调用回调的同时，才会结束缓冲操作
 */
const confirm = async (newV, callback) => {
  emits('modify', newV, () => {
    showModal.value = false
    callback()
  })
}
</script>

<style lang="scss" scoped></style>
