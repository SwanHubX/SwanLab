<template>
  <div class="w-full h-full py-5">
    <!-- 导航栏 -->
    <div class="px-6 border-b relative">
      <!-- 第一行内容，项目标题、实验标题、编辑按钮、删除按钮 -->
      <div
        class="flex items-center gap-3 transition-transform duration-300"
        :class="{ 'translate-x-8': !isSideBarShow }"
      >
        <!-- 项目标题/实验标题 -->
        <h1 class="text-2xl">
          <RouterLink class="hover:underline underline-offset-2" to="/">{{ projectStore.name }}</RouterLink>
          /
          <span class="font-semibold">{{ experimentStore.name }}</span>
        </h1>
        <ConfigEditor type="experiment" @modify="modifyExperiment" :disabled="experimentStore.isRunning" />
        <DeleteButton
          class="absolute right-6"
          type="experiment"
          :disabled="experimentStore.isRunning"
          @confirm="deleteExperiment"
        />
      </div>
      <!-- 第二行内容，实验描述 -->
      <p class="experiment-description" v-if="experimentStore.description">
        {{ experimentStore.description }}
      </p>
      <!-- 第三行内容，导航标签 -->
      <nav class="flex items-center gap-8">
        <RouterLink
          class="nav-item"
          active-class="nav-active"
          :data-text="nav.label"
          v-for="nav in navs"
          :key="nav.to"
          :to="nav.to"
        >
          {{ nav.label }}
        </RouterLink>
      </nav>
    </div>
    <div class="experiment-content">
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
const router = useRouter()
const projectStore = useProjectStore()
const experimentStore = useExperimentStroe()
const showErrorView = inject('showErrorView')

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
      showErrorView(data.code, data.message)
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
    to: `/experiment/${experimentStore.id}/index`
  },
  {
    label: t('experiment.navs.chart'),
    to: `/experiment/${experimentStore.id}/chart`
  },
  {
    label: t('experiment.navs.log'),
    to: `/experiment/${experimentStore.id}/log`
  },
  {
    label: t('experiment.navs.env'),
    to: `/experiment/${experimentStore.id}/env`
  }
]
</script>

<style scoped lang="scss">
.experiment-content {
  @apply w-full overflow-y-auto;
}

.experiment-description {
  @apply mt-3.5 w-full break-words text-sm h-10;
  display: -webkit-box;
  -webkit-box-orient: vertical;
  overflow: hidden;
  -webkit-line-clamp: 2; /* 设置为希望显示的最大行数 */
}

.nav-item {
  @apply px-2.5 py-2 relative text-lg text-dimmer;
}
.nav-active {
  @apply text-positive-higher;
  &:after {
    @apply w-full h-0.5 bg-positive-higher absolute -bottom-[1px] left-1/2 -translate-x-1/2;
    content: '';
  }
  // 字体加粗
  // &::before {
  //   @apply absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2 font-semibold whitespace-nowrap;
  //   content: attr(data-text);
  //   color: var(--positive-higher);
  // }
}
</style>
