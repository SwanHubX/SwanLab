<template>
  <div class="w-full h-full py-5">
    <!-- 导航栏 -->
    <div class="px-6 border-b">
      <!-- 第一行内容，项目标题、实验标题、编辑按钮、删除按钮 -->
      <div class="experiment-title transition-marging duration-300">
        <div class="flex items-center gap-3">
          <!-- 项目标题/实验标题 -->
          <h1 class="text-2xl items-center gap-1 truncate max-w-sm sm:max-w-lg 2xl:max-w-5xl">
            <span class="font-semibold">{{ experimentStore.name }}</span>
          </h1>
          <!-- 编辑按钮 -->
          <ConfigEditor type="experiment" @modify="modifyExperiment" :disabled="experimentStore.isRunning" />
          <!-- 实验状态 -->
          <SLStatusLabel :name="experiment.name" :status="experiment.status">
            {{ $t('experiment.status.' + experiment.status) }}
          </SLStatusLabel>
          <slot name="stop-button" v-if="experimentStore.isRunning"></slot>
        </div>
        <!-- 删除按钮 -->
        <div class="flex justify-end grow transition-padding duration-300 ml-1">
          <DeleteButton type="experiment" @confirm="deleteExperiment" />
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
    <slot></slot>
  </div>
</template>

<script setup>
/**
 * @description: 实验页布局（不含左侧侧边栏）
 * @file: ExperimentLayout.vue
 * @since: 2023-12-09 20:22:32
 **/
import ConfigEditor from '@swanlab-vue/business/config-editor/ConfigEditor.vue'
import DeleteButton from '@swanlab-vue/business/config-editor/DeleteButton.vue'
import http from '@swanlab-vue/api/http'
import { useProjectStore, useExperimentStore } from '@swanlab-vue/store'
import { useRouter } from 'vue-router'
import { message } from '@swanlab-vue/components/message'
import { t } from '@swanlab-vue/i18n'
import { ref } from 'vue'

const router = useRouter()
const projectStore = useProjectStore()
const experimentStore = useExperimentStore()
const experiment = ref(experimentStore.experiment)

// ---------------------------------- 删除实验 ----------------------------------
const deleteExperiment = () => {
  http
    .delete(`/experiment/${experimentStore.id}`)
    .then(({ data }) => {
      projectStore.deleteExperiment(data.experiment_id)
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
  data.name = data.name.trim()
  experimentStore.setExperimentPartial(data)
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
    icon: 'chart'
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
  @apply mt-3 w-full break-words text-sm text-dimmer;
  display: -webkit-box;
  -webkit-box-orient: vertical;
  overflow: hidden;
  -webkit-line-clamp: 2; /* 设置为希望显示的最大行数 */
}

.experiment-navs {
  @apply flex items-center gap-8 mt-3 w-full overflow-x-auto overflow-y-visible;
  // 隐藏滚动条
  &::-webkit-scrollbar {
    display: none;
  }
  .nav-item {
    @apply px-1 pt-2 pb-1.5 relative text-dimmer whitespace-nowrap ring-0 outline-none;
    @apply mb-1 rounded flex items-center gap-1.5;
    &:hover {
      @apply bg-higher;
    }
  }
  .nav-active {
    @apply text-default font-semibold;
    &:hover {
      background-color: transparent !important;
    }
    &:after {
      @apply w-full h-0.5 bg-positive-higher absolute -bottom-1 left-1/2 -translate-x-1/2 z-10;
      content: '';
    }
  }
}
</style>
