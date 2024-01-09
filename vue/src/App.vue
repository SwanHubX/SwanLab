<template>
  <MainLayout :version="version" v-if="ready">
    <router-view v-if="!errorCode" />
  </MainLayout>
  <ErrorView :code="errorCode" :message="errorMessage" v-if="errorCode" />
  <!-- 全局气泡提示 -->
  <SLMessages ref="messagesRef" />
  <!-- 全局确认弹窗 -->
  <SLComfirm ref="confirmRef" />
</template>

<script setup>
import MainLayout from './layouts/main/MainLayout.vue'
import ErrorView from './views/error/ErrorView.vue'
import http from './api/http'
import { useProjectStore } from '@swanlab-vue/store'
import { computed } from 'vue'
import { ref } from 'vue'
import { provide } from 'vue'
import { useRoute } from 'vue-router'
import { watch } from 'vue'
import { installMessage, SLMessages, message } from '@swanlab-vue/components/message'
import { installConfirm, SLComfirm } from './components/comfirm'
import { onMounted } from 'vue'

const projectStore = useProjectStore()
const ready = ref()

// ---------------------------------- 在此处请求项目信息 ----------------------------------
const version = ref()
http
  .get('/project')
  .then(({ data, _header }) => {
    projectStore.setProject(data)
    version.value = _header['swanlab-version']
  })
  .catch((response) => {
    // console.error(response)
    errorCode.value = response.data?.code || 3000 // 3000 时，后端启动失败
    version.value = response.headers['swanlab-version']
  })
  .finally(() => {
    ready.value = true
  })

// ---------------------------------- 错误处理 ----------------------------------

const errorCode = ref(0) // 错误码
const errorMessage = ref('') // 错误信息

/**
 * 跳转错误页
 * @param {number} code 错误码
 * @param {string} message 错误信息，可选，在没有自定义错误页的时候使用
 */
const showErrorView = (code, message = '') => {
  errorCode.value = code
  errorMessage.value = message
}

/**
 * 清除错误缓存
 */
const clearErrorView = () => {
  errorCode.value = 0
  errorMessage.value = ''
}

provide('showErrorView', showErrorView)
// provide('clearErrorView', clearErrorView)

const route = useRoute()

// 监测路由修改
watch(
  computed(() => route.fullPath),
  () => {
    // 清除错误，恢复正常页面
    clearErrorView()
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
