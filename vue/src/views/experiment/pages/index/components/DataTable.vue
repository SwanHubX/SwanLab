<template>
  <div class="table">
    <!-- title and search -->
    <div class="flex items-center justify-between border-b py-2 px-4 bg-higher">
      <p class="pr-10">{{ title }}</p>
      <SLSearch
        class="bg-default max-w-52"
        :placeholder="$t('experiment.index.config.table.search')"
        reverse
        v-model="search"
      />
    </div>
    <!-- table -->
    <div class="w-full overflow-x-auto">
      <!-- header -->
      <div class="line bg-higher">
        <div class="header-item" v-for="item in column" :key="item.key" :title="item.title">{{ item.title }}</div>
      </div>
      <!-- body -->
      <TableLine class="line" v-for="line in tableData" :key="line.key" :line="line">
        <div class="body-item hover:before:contents" v-for="item in line" :key="item.key" :title="item">
          {{ item }}
        </div>
      </TableLine>
      <!-- empty body -->
      <div class="w-full py-4 flex justify-center" v-if="data.length === 0">
        {{ $t('experiment.index.config.table.empty') }}
      </div>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 展示实验配置、实验总结数据的表格
 * @file: DataTable.vue
 * @since: 2024-01-24 14:40:21
 **/

import { ref, computed } from 'vue'
import TableLine from './TableLine.vue'

// ---------------------------------- 组件接口 ----------------------------------

const props = defineProps({
  title: {
    type: String,
    default: ''
  },
  // 表格体的数据
  data: {
    type: Array,
    default: () => {
      return []
    }
  },
  // 表头
  column: {
    type: Array,
    default: () => {
      return []
    }
  }
})

// ---------------------------------- 过滤 ----------------------------------

const search = ref('')

const tableData = computed(() => {
  if (search.value === '') {
    return props.data
  }
  return props.data.filter((item) => {
    return Object.keys(item).some((key) => {
      return String(item[key]).toLowerCase().includes(search.value.toLowerCase())
    })
  })
})
</script>

<style lang="scss" scoped>
.table {
  @apply w-full border rounded-lg overflow-hidden;
}

.line {
  @apply grid grid-cols-2 gap-3;
}

.header-item {
  @apply w-full border-r pl-4 my-3;
  &:last-child {
    @apply border-none pl-0;
  }
}

.body-item {
  @apply py-3 pl-4 pr-2 truncate;
  &:last-child {
    @apply border-none pl-0;
  }
}
</style>
