<template>
  <div class="w-full h-full py-5 relative">
    <!-- 导航栏 -->
    <div class="px-6 border-b">
      <!-- 第一行内容，项目标题、实验标题、编辑按钮、删除按钮 -->
      <div class="experiment-title transition-marging duration-300" :class="{ 'ml-8': !isSideBarShow }">
        <div class="flex items-center gap-3">
          <!-- 项目标题/实验标题 -->
          <h1 class="text-2xl items-center gap-1 truncate max-w-sm sm:max-w-lg 2xl:max-w-5xl">
            <RouterLink class="hover:underline underline-offset-2" to="/">{{ projectStore.name }}</RouterLink>
            /
            <span class="font-semibold">{{ experimentStore.name }}</span>
          </h1>
          <!-- 编辑按钮 -->
          <ConfigEditor type="experiment" @modify="modifyExperiment" :disabled="experimentStore.isRunning" />
          <!-- 实验状态 -->
          <SLStatusLabel :name="experiment.name" :id="experiment.id" :status="experiment.status" />
          <slot name="stop-button" v-if="experimentStore.isRunning"></slot>
        </div>
        <!-- 删除按钮 -->
        <div class="flex justify-end grow transition-padding duration-300 ml-1" :class="{ 'pr-8': !isSideBarShow }">
          <DeleteButton type="experiment" :disabled="experimentStore.isRunning" @confirm="deleteExperiment" />
        </div>
      </div>
      <!-- 第二行内容，实验描述 -->
      <p class="experiment-description" v-if="experimentStore.description">
        {{ experimentStore.description }}
      </p>
      <!-- 第三行内容，导航标签 -->
      <nav class="experiment-navs">
        <RouterLink
          class="nav-item"
          active-class="nav-active"
          :data-text="nav.label"
          v-for="nav in navs"
          :key="nav.to"
          :to="nav.to"
        >
          <SLIcon :icon="nav.icon" class="w-5 h-5" />
          {{ nav.label }}
        </RouterLink>
      </nav>
    </div>
    <div class="w-full overflow-y-auto pb-10">
      <slot></slot>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 实验页布局（不含左侧侧边栏）
 * @file: ExperimentLayout.vue
 * @since: 2023-12-09 20:22:32
 **/
import ConfigEditor from '@swanlab-vue/components/config-editor/ConfigEditor.vue'
import DeleteButton from '@swanlab-vue/components/config-editor/DeleteButton.vue'
import http from '@swanlab-vue/api/http'
import { useProjectStore, useExperimentStroe } from '@swanlab-vue/store'
import { inject } from 'vue'
import { useRouter } from 'vue-router'
import { message } from '@swanlab-vue/components/message'
import { t } from '@swanlab-vue/i18n'
import { ref } from 'vue'

const router = useRouter()
const projectStore = useProjectStore()
const experimentStore = useExperimentStroe()
const experiment = ref(experimentStore.experiment)

// ---------------------------------- 控制h1缩进 ----------------------------------
const isSideBarShow = inject('isSideBarShow')

// ---------------------------------- 删除实验 ----------------------------------
const deleteExperiment = () => {
  http
    .delete(`/experiment/${experimentStore.id}`)
    .then(({ data }) => {
      projectStore.setProject(data.project)
      router.replace('/').then(() => {
        message.success('Delete Successfully')
      })
    })
    .catch(({ data }) => {
      message.error(data.message)
    })
}

// ---------------------------------- 修改实验信息 ----------------------------------

const modifyExperiment = async (newV, hideModal) => {
  const id = experimentStore.id
  const { data } = await http.patch(`/experiment/${id}`, newV)
  experimentStore.setExperiment(data.experiment)
  projectStore.setExperimentInfo(id, newV)
  hideModal()
}

// ---------------------------------- 导航标签配置 ----------------------------------
const navs = [
  {
    label: t('experiment.navs.index'),
    to: `/experiment/${experimentStore.id}/index`,
    icon: 'experiment'
  },
  {
    label: t('experiment.navs.chart'),
    to: `/experiment/${experimentStore.id}/chart`,
    icon: 'charts'
  },
  {
    label: t('experiment.navs.log'),
    to: `/experiment/${experimentStore.id}/log`,
    icon: 'logs'
  },
  {
    label: t('experiment.navs.env'),
    to: `/experiment/${experimentStore.id}/env`,
    icon: 'info'
  }
]
</script>

<style scoped lang="scss">
.experiment-title {
  @apply flex items-center w-full overflow-x-auto overflow-y-visible;
  // 隐藏滚动条
  &::-webkit-scrollbar {
    display: none;
  }
}
.experiment-description {
  @apply mt-3.5 w-full break-words text-sm text-dimmer;
  display: -webkit-box;
  -webkit-box-orient: vertical;
  overflow: hidden;
  -webkit-line-clamp: 2; /* 设置为希望显示的最大行数 */
}

.experiment-navs {
  @apply flex items-center gap-8 mt-6 w-full overflow-x-auto overflow-y-visible -ml-3;
  // 隐藏滚动条
  &::-webkit-scrollbar {
    display: none;
  }
  .nav-item {
    @apply px-2.5 pt-2 pb-1.5 relative text-lg text-dimmer whitespace-nowrap ring-0 outline-none;
    @apply mb-1 rounded flex items-center gap-1.5;
    &:hover {
      @apply bg-higher;
    }
  }
  .nav-active {
    @apply text-positive-higher;
    &:hover {
      background-color: transparent !important;
    }
    &:after {
      @apply w-full h-0.5 bg-positive-higher absolute -bottom-1 left-1/2 -translate-x-1/2 z-10;
      content: '';
    }
    // 字体加粗
    // &::before {
    //   @apply absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2 font-semibold whitespace-nowrap;
    //   content: attr(data-text);
    //   color: var(--positive-higher);
    // }
  }
}
</style>
