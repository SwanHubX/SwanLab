<template>
  <div class="w-full h-full bg-higher flex justify-between items-center px-6">
    <!-- logo and version -->
    <div class="flex">
      <SLIcon class="w-8 h-8" icon="logo" />
      <div class="flex items-end pb-0.5">
        <span class="font-semibold text-2xl leading-6 mr-0.5">Swanlab</span>
        <span class="whitespace-nowrap text-xs pl-2 text-dimmer"> v{{ version }}</span>
      </div>
    </div>
    <!-- links -->
    <div class="w-full h-full px-8 flex items-center font-semibold">
      <a class="hover:pb-1 transition-all" :href="item.link" target="_blank" v-for="item in links" :key="item.link">{{
        item.title
      }}</a>
    </div>
    <!-- fixeds -->
    <div class="flex items-center font-semibold">
      <a
        :href="item.link"
        target="_blank"
        class="flex gap-1.5 items-center h-full"
        v-for="item in fixeds"
        :key="item.icon"
      >
        <SLIcon :icon="item.icon" class="h6 w-6" />
        {{ item.title }}
      </a>
    </div>
  </div>
</template>

<script setup>
/**
 * @description:
 * @file: PageHeader.vue
 * @since: 2024-01-09 11:13:20
 **/

import { ref } from 'vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'

defineProps({
  version: {
    type: String,
    default: 'unknown'
  }
})

// ---------------------------------- 链接配置 ----------------------------------

const links = ref([
  {
    title: 'Docs',
    link: 'https://geektechstudio.feishu.cn/wiki/space/7310593325374013444?ccm_open_type=lark_wiki_spaceLink&open_tab_from=wiki_home'
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

<style lang="scss" scoped></style>
