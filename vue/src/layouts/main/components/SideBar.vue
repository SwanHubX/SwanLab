<template>
  <!-- 剩余区域 -->
  <div class="flex flex-col grow h-full">
    <!-- 概览区域 -->
    <div class="p-4 flex flex-col border-b pr-11">
      <RouterLink to="/" active-class="active-router">
        <SLIcon icon="dashboard" class="w-4 h-4 mr-2" />
        <!-- <span>{{ $t('sider.nav.home') }}</span> -->
        <span>Project Dashboard</span>
      </RouterLink>
    </div>
    <!-- 实验路由 -->
    <div class="experiments-container">
      <SLSearch @input="search" />
      <!-- 实验列表 -->
      <RouterLink
        v-for="experiment in experiments"
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
  </div>
</template>

<script setup>
/**
 * @description: 侧边栏导航
 * @file: HomeSiderBar.vue
 * @since: 2023-12-04 18:20:02
 **/
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import SLSearch from '@swanlab-vue/components/SLSearch.vue'
import { ref, computed } from 'vue'
import { RouterLink } from 'vue-router'
import { useProjectStore } from '@swanlab-vue/store'

const projectStore = useProjectStore()
// ---------------------------------- 实验id转路由 ----------------------------------
const getExperimentRouter = (experiment) => {
  return `/experiment/${experiment.experiment_id}`
}

// ---------------------------------- 搜索实验 ----------------------------------

// 需要展示的实验信息——默认展示全部，但在搜索过后，更新为搜索结果
const experiments = computed(() => {
  if (!searchValue.value) return projectStore.experiments
  return projectStore.experiments.filter((expr) => expr.name.toLowerCase().includes(searchValue.value))
})

const searchValue = ref('')

const search = (value) => {
  searchValue.value = value.toLowerCase()
}
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
