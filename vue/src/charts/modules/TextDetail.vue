<template>
  <div class="px-6">
    <!-- infos -->
    <div v-for="item in info_list" :key="item" class="flex items-center pb-3">
      <span class="block w-48 font-semibold">{{ $t(`common.chart.text-chart.titles.${item.key}`) }}:</span>
      <span>{{ item.value }}</span>
    </div>
    <!-- text -->
    <p class="font-semibold pb-2">{{ $t('common.chart.text-chart.titles.text') }}</p>
    <div class="p-4 min-h-[200px] border rounded bg-higher">{{ data.line.data }}</div>
  </div>
</template>

<script setup>
/**
 * @description: 展示 text 详细信息，与弹窗中调用
 * @file: TextDetail.vue
 * @since: 2024-02-21 13:46:40
 **/
import { computed } from 'vue'

const props = defineProps({
  data: {
    type: Object,
    required: true,
    default: () => {}
  }
})

const info_list = computed(() => {
  return [
    {
      key: 'tag',
      value: props.data.tag
    },
    {
      key: 'step',
      value: props.data.line.index
    },
    {
      key: 'caption',
      value: props.data.line.more.caption
    },
    {
      key: 'count',
      value: props.data.line.data.length
    }
  ]
})
</script>

<style lang="scss" scoped></style>
