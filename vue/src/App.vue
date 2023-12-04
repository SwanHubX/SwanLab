<template>
  <MainLayout>
    <template #left>
      <HomeSiderBar />
    </template>
    <router-view />
  </MainLayout>
</template>

<script setup>
import MainLayout from './layouts/MainLayout.vue'
import HomeSiderBar from './components/HomeSiderBar.vue'
</script>

<style scoped></style>
