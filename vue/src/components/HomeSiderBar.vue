<template>
  <div class="w-full h-full">
    <!-- 头部卡片 -->
    <div class="pt-4 pb-2 flex justify-between items-center px-4">
      <!-- 项目信息 -->
      <div class="flex items-end">
        <SlabIcon icon="logo" class="w-6 h-6 mr-1" />
        <span class="font-bold text-2xl leading-6 pr-0.5">Swanlab</span>
      </div>
      <!-- 友链 -->
      <div class="flex">
        <a href="https://github" class="link">
          <SlabIcon icon="github" class="w-5 h-5" />
        </a>
        <a href="https://swanhub.co" class="link">
          <SlabIcon icon="logo" class="w-5 h-5" />
        </a>
      </div>
    </div>
    <!-- 版本 -->
    <div class="px-4 pb-4 text-xs text-dimmer border-b border-default">
      <span>{{ $t('sider.version.version', { version: '1.0.0' }) }}</span>
      <a href="https://swanhub.co" class="hover:underline">{{ $t('sider.version.updates') }}</a>
    </div>
    <!-- 概览 -->
    <RouterLink to="/" class="router-link my-4 mx-4" active-class="active-router">
      <SlabIcon icon="home" class="w-4 h-4 mr-2 text-positive-default" />
      <span>{{ $t('sider.nav.home') }}</span>
    </RouterLink>
    <!-- 实验路由 -->
    <div class="px-4 flex flex-col pt-4 border-t border-default">
      <input type="text" />
      <RouterLink
        v-for="experiment in display_experiments"
        :key="experiment.experiment_id"
        :to="`/experiment/${experiment.experiment_id}`"
        class="router-link"
        activeClass="active-router"
      >
        <SlabIcon icon="experiment" class="w-4 h-4 mr-3" />
        <span>{{ experiment.name }}</span>
      </RouterLink>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 侧边栏导航
 * @file: HomeSiderBar.vue
 * @since: 2023-12-04 18:20:02
 **/
import SlabIcon from './SlabIcon.vue'
import { ref } from 'vue'
import http from '@swanlab-vue/api/http'

// ---------------------------------- 获取已有实验列表 ----------------------------------
// 所有的实验信息
const experiments = ref([])
// 需要展示的实验信息——默认展示全部，但在搜索过后，更新为搜索结果
const display_experiments = ref([])
http
  .get('/mock/api/project/experiments')
  .then((res) => {
    console.log(res.data)
    experiments.value = res.data.experiments
    display_experiments.value = res.data.experiments
  })
  .catch((error) => {
    console.error(error)
  })

// ---------------------------------- 搜索实验 ----------------------------------
</script>

<style lang="scss" scoped>
.link {
  @apply p-1.5 bg-default ml-2 rounded border border-default transition-all hover:rounded-full;
}

.router-link {
  @apply flex items-center px-4 py-3 text-default;
}

.active-router {
  @apply bg-positive-highest text-positive-default font-semibold rounded-lg;
}
</style>
