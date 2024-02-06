<template>项目对比图表 TODO...</template>

<script setup>
import http from '@swanlab-vue/api/http'

/**
 * @description: 项目对比图表，本组件完成项目对比图表的数据的请求和展示
 * 数据请求分为
 * @file: ChartsView.vue
 * @since: 2024-01-27 13:05:27
 **/

http.get('/project/charts').then((res) => {
  console.log(res)
})
</script>

<style lang="scss" scoped></style>
