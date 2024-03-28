<template>
  <Dialog as="div" class="sl-dialog-mask" data-sl-modal :open="modelValue" @click.prevent="close('mask')">
    <!-- 模态框主体部分 -->
    <div class="sl-dialog animate-bounce-in" :class="$props.class" :style="{ maxWidth }" data-sl-modal @click.stop="">
      <!-- 关闭按钮 -->
      <button
        class="absolute z-full !border-none"
        :class="[buttonPosition, buttonStyle]"
        @click="close('button')"
        v-if="!hiddenClose"
      >
        <SLIcon class="w-full h-full" icon="close" />
      </button>
      <!-- 插槽，模态框内容 -->
      <slot></slot>
    </div>
  </Dialog>
</template>

<script setup>
import { computed, watch } from 'vue'
import { Dialog } from '@headlessui/vue'
import SLIcon from './SLIcon.vue'

const props = defineProps({
  // 组件挂载时是否显示，默认不显示，值变为true时显示
  modelValue: {
    type: Boolean,
    default: false
  },
  // 点击遮罩层是否关闭，默认不关闭
  closeOnOverlayClick: {
    type: Boolean,
    default: false
  },
  // 关闭按钮的位置
  buttonPosition: {
    type: String,
    default: 'right-4 top-4 '
  },
  // 除了位置以外的其他样式
  buttonStyle: {
    type: String,
    default: 'w-5 h-5 p-1 rounded'
  },
  maxW: {
    type: String,
    default: '800'
  },
  class: {
    type: String,
    default: ''
  },
  hiddenClose: {
    type: Boolean,
    default: false
  },
  // 是否可通过键盘 esc 键关闭弹窗
  escExit: {
    type: Boolean,
    default: false
  }
})
// onBeforeClose事件特指人为关闭事件，如果是在外部用js修改modelValue的值，不会触发该事件
const emits = defineEmits(['onBeforeClose', 'update:modelValue', 'onExit'])

// 关闭弹窗方法
const close = (type) => {
  if (type === 'mask' && !props.closeOnOverlayClick) return
  emits('update:modelValue', false)
  emits('onBeforeClose')
}

const maxWidth = computed(() => {
  return props.maxW + 'px'
})

// ---------------------------------- 键盘esc关闭弹窗 ----------------------------------

/**
 * 通过按键关闭弹窗
 * @param {*} param0
 */
const handleEsc = ({ key }) => {
  if (key != 'Escape') return
  emits('onExit')
  if (props.escExit) {
    emits('onBeforeClose')
    close('esc')
  }
}

// 订阅或销毁键盘事件
watch(
  () => props.modelValue,
  (value) => {
    if (value) return window.addEventListener('keyup', handleEsc)
    window.removeEventListener('keyup', handleEsc)
  }
)
</script>

<style lang="scss">
.sl-dialog-mask[data-sl-modal] {
  @apply z-full w-full h-full fixed top-0 left-0;
  @apply py-16 px-8 overflow-y-auto overflow-x-hidden;
  background-color: var(--background-overlay);
  // 高斯模糊
  backdrop-filter: blur(10px);
}

.sl-dialog[data-sl-modal] {
  @apply rounded-lg relative mx-auto;
  /** 样式 */
  background-color: var(--background-default);
  border: 1px solid var(--outline-default);
}
</style>
