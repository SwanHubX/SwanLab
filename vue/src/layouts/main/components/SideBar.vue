<template>
  <!-- 侧边栏区域 -->
  <div class="flex flex-col grow h-full bg-higher">
    <!-- 概览区域 -->
    <div class="p-4 flex flex-col border-b gap-2 font-semibold">
      <!-- 项目信息 -->
      <h1 class="font-semibold mb-1">{{ $t('common.sidebar.project.title') }}</h1>
      <RouterLink to="/" active-class="active-link">
        <SLIcon icon="runs" class="w-4 h-4 mr-2" />
        <span>{{ $t('common.sidebar.project.runs') }}</span>
      </RouterLink>
      <RouterLink to="/charts" active-class="active-link">
        <SLIcon icon="charts" class="w-4 h-4 mr-2" />
        <span>{{ $t('common.sidebar.project.charts') }}</span>
      </RouterLink>
    </div>
    <!-- 实验路由 -->
    <div class="experiments-container" ref="expContainerRef">
      <div class="flex flex-col gap-4 sticky top-0 bg-higher py-4 -mb-4">
        <div class="flex items-center gap-1.5">
          <h1 class="font-semibold">{{ $t('common.sidebar.experiments.title') }}</h1>
          <span class="font-semibold px-3 text-sm bg-highest rounded-full">{{ totalExperiments }}</span>
        </div>
        <SLSearch @input="search" reverse :placeholder="$t('common.search.placeholder')" />
      </div>
      <!-- 实验列表 -->
      <RouterLink
        v-for="experiment in experiments"
        :key="experiment.id"
        :to="getExperimentRouter(experiment)"
        :title="experiment.name"
        class="experiment-link"
        active-class="active-link"
      >
        <circle
          class="w-4 h-4 rounded-full mr-3 flex-shrink-0"
          :style="{ backgroundColor: getExperimentColor(experiment) }"
        />
        <span class="truncate">{{ experiment.name }}</span>
        <!-- 更多信息，进入此容器后不触发父容器所有效果(包括hover、active的css效果) -->
        <div class="more-info" @click.prevent @mouseenter="removeHover" @mouseleave="resetHover">
          <!-- 如果实验正在运行，显示running -->
          <span v-if="experiment.status === 0"> ({{ $t('common.sidebar.experiments.running') }}) </span>
          <!-- 如果在charts页面，显示眼睛 -->
          <button class="show-button" v-if="$route.path == '/charts'" @click="changeExperimentShow(experiment.id)">
            <SLIcon icon="eye" class="w-full h-full" v-if="experiment.show" />
            <SLIcon icon="eye-close" class="w-full h-full text-dimmest" v-else />
          </button>
        </div>
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
import http from '@swanlab-vue/api/http'
import { debounces } from '@swanlab-vue/utils/common'

const projectStore = useProjectStore()
// ---------------------------------- 实验id转路由 ----------------------------------
const getExperimentRouter = (experiment) => {
  return `/experiment/${experiment.id}`
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

// ---------------------------------- 获取实验颜色 ----------------------------------
const getExperimentColor = (experiment) => {
  return experiment.light
}

// ---------------------------------- 计算实验数量，也包括可视实验数量 ----------------------------------

const totalExperiments = computed(() => {
  return projectStore.experiments.length
})

// ---------------------------------- 项目图表界面下，点击眼睛后的效果 ----------------------------------
const path = debounces((id, show) => {
  http.patch('/experiment/' + id + '/show', {
    show
  })
}, 500)
const changeExperimentShow = (id) => {
  const show = projectStore.changeExperimentShow(id)
  // 请求后端更新，直接更新store后请求，不需要等待
  path(id, show)
}

// ---------------------------------- 进入more-info部分，将父元素hover效果移除 ----------------------------------
const removeHover = (e) => {
  // console.log('进入', e.target.parentNode)
  const parent = e.target.parentNode
  parent.classList.add('!bg-transparent')
}
const resetHover = (e) => {
  // console.log('离开', e.target.parentNode)
  const parent = e.target.parentNode
  parent.classList.remove('!bg-transparent')
}
</script>

<style lang="scss" scoped>
// RouterLink共享样式
a {
  @apply flex items-center px-4 h-9 text-default hover:bg-highest rounded-lg;
}

.active-link {
  @apply bg-highest text-default cursor-default;
}

.experiments-container {
  @apply flex flex-col p-4 pt-0 grow gap-4 overflow-auto;
  &::-webkit-scrollbar-track {
    background: transparent;
  }

  // 实验链接样式
  .experiment-link {
    @apply flex-shrink-0 text-base pr-0;
    .more-info {
      @apply flex justify-end grow text-primary-default pr-4 gap-2 cursor-default;
      .show-button {
        @apply w-7 h-6 p-1 rounded;
        &:hover {
          @apply text-primary-highest bg-highest;
        }
      }
    }
  }
}
</style>
