<template>
  <div class="w-full h-full bg-dimmest text-dimmest flex justify-between items-center px-6">
    <!-- logo and version -->
    <div class="flex items-center gap-1.5 hover:cursor-pointer" @click="$router.push('/')">
      <!-- icon -->
      <HeaderIcon />
      <!-- version -->
      <div class="flex items-end">
        <span class="font-semibold mr-0.5">SwanLab</span>
        <!-- 版本号被注释了 -->
        <!-- <span class="whitespace-nowrap text-xs pl-2 text-dimmer"> {{ formatVersion(version) }}</span> -->
      </div>
    </div>
    <!-- links -->
    <div class="w-full grow px-8 items-center font-semibold justify-end gap-8 md:flex hidden">
      <a :href="item.link" target="_blank" class="hover:text-white-higher" v-for="item in links" :key="item.link">
        {{ item.title }}
      </a>
    </div>
    <!-- fixeds -->
    <div class="flex items-center font-semibold">
      <a
        :href="item.link"
        target="_blank"
        class="flex gap-1.5 items-center h-full text-dimmest hover:text-white-higher"
        v-for="item in fixeds"
        :key="item.icon"
      >
        <SLIcon :icon="item.icon" class="h-8 w-8" />
        <!-- {{ item.title }} -->
      </a>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 顶部页头
 * @file: PageHeader.vue
 * @since: 2024-01-09 11:13:20
 **/

import { ref } from 'vue'
import HeaderIcon from './HeaderIcon.vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'

defineProps({
  version: {
    type: String,
    default: 'unknown'
  }
})

// ---------------------------------- 格式化版本号 ----------------------------------
// const formatVersion = (version) => {
//   if (version === 'unknown') return version
//   return 'v' + version
// }

// ---------------------------------- 链接配置 ----------------------------------

const links = ref([
  {
    title: 'Docs',
    link: 'https://geektechstudio.feishu.cn/wiki/space/7310593325374013444?ccm_open_type=lark_wiki_spaceLink&open_tab_from=wiki_home'
  },
  {
    title: 'Examples',
    link: 'https://github.com/SwanHubX/SwanLab-examples'
  },
  {
    title: 'Feedback',
    link: 'https://github.com/SwanHubX/SwanLab/issues'
  }
])

const fixeds = ref([
  {
    title: 'GitHub',
    icon: 'github',
    link: 'https://github.com/SwanHubX/SwanLab'
  }
])
</script>

<style lang="scss" scoped>
.a-hover {
  &:hover {
    @apply text-white-higher;
  }
}
</style>
