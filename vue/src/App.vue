<template>
  <MainLayout v-if="ready">
    <template #left>
      <HomeSiderBar />
    </template>
    <router-view v-if="!error_code" />
    <ErrorView :code="error_code" v-else />
  </MainLayout>
</template>

<script setup>
import MainLayout from './layouts/MainLayout.vue'
import HomeSiderBar from './components/HomeSiderBar.vue'
import ErrorView from './views/error/ErrorView.vue'
import http from './api/http'
import { useProjectStore } from '@swanlab-vue/store'
import { computed } from 'vue'
import { ref } from 'vue'
import { provide } from 'vue'
import { useRoute } from 'vue-router'
import { watch } from 'vue'
const projectStore = useProjectStore()
const ready = computed(() => !!projectStore.experiments || error_code.value)

// ---------------------------------- 在此处请求项目信息 ----------------------------------
http
  .get('/project')
  .then(({ data }) => {
    projectStore.setProject(data)
  })
  .catch((response) => {
    // console.error(response)
    error_code.value = response.data?.code || 3000 // 3000 时，后端启动失败
  })

// ---------------------------------- 错误处理 ----------------------------------
const error_code = ref(0)
const show_error = (code) => {
  error_code.value = code
}
const clear_error = () => {
  error_code.value = 0
}
provide('show_error', show_error)
provide('clear_error', clear_error)

const route = useRoute()
watch(
  computed(() => route.fullPath),
  () => clear_error()
)
</script>

<style scoped></style>
