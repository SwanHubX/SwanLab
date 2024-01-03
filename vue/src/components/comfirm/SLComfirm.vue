<template>
  <SLModal @on-before-close="cancel" v-model="visiable" class="p-6 h-48 flex flex-col" button-position="top-7 right-6">
    <h1 class="text-xl font-semibold">{{ text.title }}</h1>
    <p class="mt-4 text-dimmer text-base">{{ text.content }}</p>
    <div class="mt-auto text-sm flex gap-3 justify-end">
      <SLButton :text="text.cancel" @click="cancel" />
      <SLButton :text="text.confirm" theme="negative" @click="confirm" />
    </div>
  </SLModal>
</template>

<script setup>
/**
 * @description: 全局确认弹窗，在app创建的、同时实例化此组件
 * @file: SLComfirm.vue
 * @since: 2024-01-03 19:19:46
 **/
import { ref } from 'vue'
import SLModal from '../SLModal.vue'
import { reactive } from 'vue'
import SLButton from '../SLButton.vue'

const visiable = ref(false)
const defaultConfig = {
  title: 'Are you sure?',
  content: 'Are you sure you want to do this?',
  confirm: 'Yes, I confirm',
  cancel: 'Cancel'
}

const text = reactive({
  title: defaultConfig.title,
  content: defaultConfig.content,
  confirm: 'Yes, I confirm',
  cancel: 'Cancel'
})

const callback = {
  resolve: () => {},
  reject: () => {}
}

// ---------------------------------- 控制展示和隐藏 ----------------------------------
// 通过调用此方法，显示弹窗，传入标题、内容、配置，配置包含：确认按钮文字、取消按钮文字、确认回调、取消回调
const show = (title, content, config, resolve, reject) => {
  visiable.value = true
  text.title = title || defaultConfig.title
  text.content = content || defaultConfig.content
  text.confirm = config?.buttonText?.confirm || defaultConfig.confirm
  text.cancel = config?.buttonText?.cancel || defaultConfig.cancel
  callback.resolve = resolve
  callback.reject = reject
}

const confirm = () => {
  visiable.value = false
  callback.resolve && callback.resolve()
}

const cancel = () => {
  visiable.value = false
  callback.reject && callback.reject()
}

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  show,
  confirm,
  cancel
})
</script>

<style lang="scss" scoped>
button {
  @apply p-2 rounded;
}
</style>
