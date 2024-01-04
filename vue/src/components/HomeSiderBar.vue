<template>
  <!-- 头部信息卡片 -->
  <div class="flex justify-between items-center p-4 border-b h-16">
    <!-- 项目信息 -->
    <div class="flex items-end">
      <SLIcon icon="logo" class="w-6 h-6 mr-1 mb-0.5" />
      <span class="font-semibold text-2xl leading-6 mr-0.5">Swanlab</span>
      <span class="whitespace-nowrap text-xs">{{ version }}</span>
    </div>
    <!-- 友链 -->
    <!-- <div class="flex gap-2">
      <a href="https://github.com" class="link" target="_blank">
        <SLIcon icon="github" class="w-4 h-4" />
      </a>
      <a href="https://swanhub.co" class="link" target="_blank">
        <SLIcon icon="logo" class="w-4 h-4" />
      </a>
    </div> -->
  </div>
  <!-- 剩余区域 -->
  <div class="flex flex-col grow max-h-[calc(100%-4rem)]">
    <!-- 概览区域 -->
    <div class="p-4 flex flex-col border-b">
      <RouterLink to="/" active-class="active-router">
        <SLIcon icon="dashboard" class="w-4 h-4 mr-2" />
        <!-- <span>{{ $t('sider.nav.home') }}</span> -->
        <span>Project Dashboard</span>
      </RouterLink>
    </div>
    <!-- 实验路由 -->
    <div class="experiments-container">
      <!-- 实验列表 -->
      <RouterLink
        v-for="experiment in projectStore.experiments"
        :key="experiment.experiment_id"
        :to="getExperimentRouter(experiment)"
        :title="experiment.name"
        class="flex-shrink-0"
        active-class="active-router"
      >
        <div class="w-4 h-4 rounded-full mr-3" :style="{ backgroundColor: experiment.color }"></div>
        <span class="truncate">{{ experiment.name }}</span>
      </RouterLink>
    </div>
    <div class="border-t border-default h-20">
      <!-- <RouterLink to="/help" class="my-4 mx-4" active-class="active-router">
        <SLIcon icon="help" class="w-4 h-4 mr-2" />
        <span>{{ $t('sider.nav.help') }}</span>
      </RouterLink> -->
      <a :href="help_url" class="my-4 mx-4" target="_blank">
        <SLIcon icon="help" class="w-4 h-4 mr-2" />
        <span>{{ $t('sider.nav.help') }}</span>
      </a>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 侧边栏导航
 * @file: HomeSiderBar.vue
 * @since: 2023-12-04 18:20:02
 **/
import SLIcon from './SLIcon.vue'
import { ref } from 'vue'
import { RouterLink } from 'vue-router'
import { useProjectStore } from '@swanlab-vue/store'

defineProps({
  version: {
    type: String,
    default: 'unknown'
  }
})

const help_url =
  'https://geektechstudio.feishu.cn/wiki/space/7310593325374013444?ccm_open_type=lark_wiki_spaceLink&open_tab_from=wiki_home'

const projectStore = useProjectStore()
// ---------------------------------- 实验id转路由 ----------------------------------
const getExperimentRouter = (experiment) => {
  return `/experiment/${experiment.experiment_id}`
}

// ---------------------------------- 搜索实验 ----------------------------------
// 需要展示的实验信息——默认展示全部，但在搜索过后，更新为搜索结果
const query = ref([])
// TODO 暂时没有实现
</script>

<style lang="scss" scoped>
.link {
  @apply p-1.5 bg-default rounded border border-default h-7;
  &:hover {
    background: var(--background-default) !important;
  }
}

a {
  @apply flex items-center px-4 h-11 text-sm text-default hover:bg-positive-dimmest rounded-lg;
}

.active-router {
  @apply bg-positive-dimmest text-positive-higher;
}

.experiments-container {
  @apply flex flex-col p-4 grow gap-2 overflow-auto;
  &::-webkit-scrollbar-track {
    background: transparent;
  }
}
</style>
