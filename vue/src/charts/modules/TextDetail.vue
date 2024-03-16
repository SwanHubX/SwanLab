<template>
  <div class="text-detail">
    <!-- infos -->
    <div v-for="item in info_list" :key="item" class="sm:flex pb-3">
      <span class="block w-48 font-semibold shrink-0">{{ $t(`chart.charts.text.titles.${item.key}`) }}:</span>
      <span :title="item.value">{{ item.value }}</span>
    </div>
    <!-- text -->
    <p class="font-semibold pb-2">{{ $t('chart.charts.text.titles.text') }}:</p>
    <div class="p-4 min-h-[20vh] max-h-[55vh] overflow-y-auto border rounded bg-default">
      <p v-for="text in data.text.split('\n')" :key="text">
        <span v-if="text != ''">{{ text }}</span>
        <br v-else />
      </p>
    </div>
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
  const line = props.data.line
  return [
    {
      key: 'tag',
      value: props.data.tag
    },
    {
      key: 'step',
      value: line.index
    },
    {
      key: 'caption',
      value: props.data.caption || '-'
    },
    {
      key: 'count',
      value: props.data.text.replace(/[\r\n\s]/g, '').length
    }
  ]
})
</script>

<style lang="scss" scoped>
.text-detail {
  @apply px-6 py-10 overflow-hidden bg-higher rounded-lg;
}
</style>
