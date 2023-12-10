<template>
  <MainLayout v-if="ready">
    <template #left>
      <HomeSiderBar />
    </template>
    <router-view />
  </MainLayout>
</template>

<script setup>
import MainLayout from './layouts/MainLayout.vue'
import HomeSiderBar from './components/HomeSiderBar.vue'
import http from './api/http'
import { useProjectStore } from '@swanlab-vue/store'
import { computed } from 'vue'
const projectStore = useProjectStore()
const ready = computed(() => !!projectStore.experiments)

// ---------------------------------- 在此处请求项目信息 ----------------------------------
;(async () => {
  const { data } = await http.get('/project')
  projectStore.setProject(data)
})()
</script>

<style scoped></style>
