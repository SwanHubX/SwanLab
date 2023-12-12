<template>
  <Tippy trigger="click">
    <template v-slot="{ hide }">
      <SLIcon icon="copy" :class="iconClass" @click="copy(hide)" />
    </template>
    <template #content> {{ $t('common.copy') }}</template>
  </Tippy>
</template>

<script setup>
/**
 * @description: 复制文字到剪切板
 * @file: SLCopy.vue
 * @since: 2023-12-12 19:46:28
 **/
import { Tippy } from 'vue-tippy'
import SLIcon from './SLIcon.vue'
import { copyTextToClipboard } from '@swanlab-vue/utils/browser'
import 'tippy.js/dist/tippy.css'

const props = defineProps({
  text: {
    // 复制到剪切板的文字
    type: String,
    required: true
  },
  iconClass: {
    type: String,
    default: 'w-4 h-4'
  },
  duration: {
    type: Number,
    default: 2000
  }
})

/**
 * 将文字复制到剪切板并显示tip
 * @param {function} hide 隐藏tip的接口
 */
const copy = (hide) => {
  copyTextToClipboard(props.text)

  // 在两秒后隐藏提示
  setTimeout(() => hide(), props.duration)
}
</script>

<style lang="scss" scoped></style>
