<template>
  <MainLayout :show-side-bar="!errorCode" v-if="ready">
    <router-view v-if="!errorCode" />
    <ErrorView :code="errorCode" :message="errorMessage" v-else />
  </MainLayout>
  <!-- 全局气泡提示 -->
  <SLMessages ref="messagesRef" />
  <!-- 全局确认弹窗 -->
  <SLConfirm ref="confirmRef" />
</template>

<script setup>
import MainLayout from './layouts/main/MainLayout.vue'
import ErrorView from './views/error/ErrorView.vue'
import http from './api/http'
import { useProjectStore } from '@swanlab-vue/store'
import { computed } from 'vue'
import { ref } from 'vue'
import { useRoute } from 'vue-router'
import { watch } from 'vue'
import { installMessage, SLMessages, message } from '@swanlab-vue/components/message'
import { installConfirm, SLConfirm } from './components/confirm'
import { onMounted } from 'vue'

// ---------------------------------- state ----------------------------------

const projectStore = useProjectStore()
const ready = ref()

// ---------------------------------- 在此处请求项目信息 ----------------------------------
http
  .get('/project')
  .then(({ data }) => {
    projectStore.setProject(data)
  })
  .catch((response) => {
    // console.error(response)
    errorCode.value = response.data?.code || 3000 // 3000 时，后端启动失败
  })
  .finally(() => {
    ready.value = true
  })

// ---------------------------------- 错误处理 ----------------------------------

const errorCode = ref(0) // 错误码
const errorMessage = ref('') // 错误信息
const route = useRoute()

// 监测路由修改
watch(
  computed(() => route.fullPath),
  (oldVal) => {
    // console.log('route change', newVal, oldVal)
    if (oldVal === undefined) return
    // 清除消息弹窗
    message.clear()
  }
)

// ---------------------------------- 项目配置 ----------------------------------

const messagesRef = ref(null)
const confirmRef = ref(null)

onMounted(() => {
  // 注册全局顶部提醒
  installMessage(messagesRef)
  installConfirm(confirmRef)
})
</script>

<style scoped></style>
