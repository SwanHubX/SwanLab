<template>
  <TransitionGroup name="list" tag="ul" class="messages-container">
    <SLMessage :message="m" v-for="m in messages" :key="m.id" />
  </TransitionGroup>
</template>

<script setup>
/**
 * @description: 消息提示组件容器，全局挂载，通过App.vue中的message对象调用
 * 此组件暴露一些方法，用于控制消息提示的显示，虽然这样破坏了vue的单向数据流，但是这样做可以减少很多重复代码
 * 并且全局理论上只会有一个SLMessages组件，所以这样做也不会有什么问题
 * @file: SLMessages.vue
 * @since: 2024-01-03 15:34:01
 **/
import { ref } from 'vue'
import { debounce } from '@swanlab-vue/utils/common'
import { uuid } from '@swanlab-vue/utils/browser'
import SLMessage from './SLMessage.vue'

const messages = ref([])
// 消息提示的最大数量，超过这个数量，添加新消息时会删除最早的消息
const maxLength = 10

/**
 * 添加一个消息提示
 * @param { String } message 消息信息
 * @param { Number } delay 消息显示的时间，单位为毫秒
 * @param { String } type 消息类型，可选值为error、success、processing
 * @param { Function } callback 消息被关闭时的回调函数，可选
 */
const add = (message, type = 'error', delay = 2000, callback) => {
  const id = uuid()
  // 生成一个函数callback，用于关闭该消息提示，并且将它传入BsMessage组件
  const close = (id, callback) => () => {
    remove(id, callback)
  }
  messages.value.push({
    id,
    text: message,
    type,
    delay,
    close: close(id, callback)
  })
  limitRemoveDebounce()
}
/**
 * 清空所有消息提示
 */
const clear = () => {
  messages.value = []
}

/**
 * 关闭一个消息提示
 * @param { String } id 消息提示的id
 */
const remove = (id, callback) => {
  // 找到该消息提示的索引
  const index = messages.value.findIndex((item) => item.id === id)
  // 如果找到了，就将该消息提示从messages中删除，没找到就不管了
  if (index !== -1) {
    messages.value.splice(index, 1)
    // 如果有回调函数，就调用回调函数
    callback && callback()
  }
}
const limitRemoveDebounce = debounce(() => {
  if (messages.value.length > maxLength) {
    const m = messages.value[0]
    m.close()
  }
}, 500)

// ---------------------------------- 暴露给外部的方法 ----------------------------------

defineExpose({
  add,
  clear
})
</script>

<style lang="scss" scoped>
.messages-container {
  @apply fixed top-10 left-1/2 -translate-x-1/2 z-message;
  @apply inline-block;
  @apply flex flex-col gap-2;
}

// vue动画
.list-enter-active,
.list-leave-active {
  transition: all 0.3s ease;
}
.list-enter-from,
.list-leave-to {
  opacity: 0;
  height: 0;
  transform: translateY(-30px);
}
</style>
