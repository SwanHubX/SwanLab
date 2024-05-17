<template>
  <div class="lc-legend">
    <div
      class="lc-legend-item"
      :style="{ color: item.color }"
      v-for="item in items"
      :key="item.name"
      @mouseenter="handleMouseenter(item)"
      @mouseleave="handleMouseleave(item)"
      @click="handleClick(item)"
    >
      <RouterLink :to="experimentPrefix + item.experiment_id">
        {{ item.name }}
      </RouterLink>
    </div>
  </div>
</template>

<script setup>
import { inject } from 'vue'

/**
 * @description: 折线图图例组件
 * @file: LineChartLegend.vue
 * @since: 2024-02-27 15:02:13
 **/

defineProps({
  items: {
    type: Array,
    default: () => []
  }
})
// console.log('items', props.items)

const experimentPrefix = inject('experimentPrefix')

const emit = defineEmits(['hoverin', 'hoverout'])

// ---------------------------------- 悬浮事件 ----------------------------------
const handleMouseenter = (item) => {
  emit('hoverin', item)
}
const handleMouseleave = (item) => {
  emit('hoverout', item)
}

// ---------------------------------- 点击事件 ----------------------------------
const handleClick = (item) => {
  emit('click', item)
}
</script>

<style lang="scss" scoped>
.lc-legend {
  max-height: 28px;
  min-height: 16px;
  @apply overflow-y-auto flex flex-wrap px-3 mb-1 gap-y-0.5 gap-x-2 leading-none text-[12px];
  @apply justify-center;
  .lc-legend-item {
    @apply flex flex-shrink-0 items-center hover:brightness-75;
    &:before {
      content: '';
      display: inline-block;
      width: 8px;
      height: 2px;
      margin-top: 2px;
      margin-right: 4px;
      background-color: currentColor;
    }
  }
}
</style>
