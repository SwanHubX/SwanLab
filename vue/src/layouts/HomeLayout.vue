<template>
  <div class="w-full h-full pt-5 flex flex-col">
    <!-- 项目标题部分 -->
    <div class="px-6 border-b flex-shrink-0">
      <!-- 第一行内容，项目标题、实验标题、编辑按钮、删除按钮 -->
      <div class="project-title transition-marging duration-300">
        <div class="flex items-center gap-3">
          <!-- 项目标题/实验标题 -->
          <h1 class="text-2xl items-center gap-1 font-semibold max-w-md truncate">
            {{ projectStore.name }}
          </h1>
          <!-- 编辑按钮 -->
          <ConfigEditor type="project" @modify="modifyProject" />
        </div>
        <!-- 删除按钮 -->
        <div class="flex justify-end grow transition-padding duration-300 ml-1">
          <DeleteButton type="project" @confirm="deleteProject" :disabled="hasRunning" />
        </div>
      </div>
      <!-- 第二行内容，项目描述 -->
      <p class="project-description" v-if="projectStore.description">
        {{ projectStore.description }}
      </p>
      <!-- 其他内容，插槽 -->
      <slot name="top"></slot>
    </div>
    <!-- 下面的内容 -->
    <slot></slot>
  </div>
</template>

<script setup>
/**
 * @description: 首页布局，主要针对首页的title和内容进行布局
 * @file: HomeLayout.vue
 * @since: 2024-01-09 16:33:56
 **/
import { computed } from 'vue'
import { useProjectStore } from '@swanlab-vue/store'
import http from '@swanlab-vue/api/http'
import { message } from '@swanlab-vue/components/message'
import ConfigEditor from '@swanlab-vue/business/config-editor/ConfigEditor.vue'
import DeleteButton from '@swanlab-vue/business/config-editor/DeleteButton.vue'
const projectStore = useProjectStore()
// ---------------------------------- 修改项目信息 ----------------------------------

const modifyProject = async (newV, hideModal) => {
  const { data } = await http.patch('/project/update', newV)
  projectStore.updateInfo(data.updates)
  hideModal()
}

// ---------------------------------- 删除项目 ----------------------------------

/**
 * 是否可以删除项目
 * 如果有实验正在进行中，直接不显示删除按钮
 */
const hasRunning = computed(() => {
  for (const experiment of projectStore.experiments) {
    if (experiment.status === 0) return true
  }
  return false
})

/**
 * 删除项目
 *
 * success：显示成功，清空状态，并显示空项目错误页
 *
 * fail：在有实验正在进行的时候删除项目
 */
const deleteProject = () => {
  http
    .delete('/project')
    .then(() => {
      message.success('Delete Successfully')
      location.reload()
    })
    .catch(() => {
      message.error('Error deleting project')
    })
}
</script>

<style lang="scss" scoped>
.project-title {
  @apply flex items-center w-full overflow-x-auto;
  // 隐藏滚动条
  &::-webkit-scrollbar {
    display: none;
  }
}

.project-description {
  @apply mt-3.5 w-full break-words text-sm text-dimmer;
  display: -webkit-box;
  -webkit-box-orient: vertical;
  overflow: hidden;
  -webkit-line-clamp: 2; /* 设置为希望显示的最大行数 */
}
</style>
