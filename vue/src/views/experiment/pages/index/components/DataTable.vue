<template>
  <div class="table">
    <!-- title and search -->
    <div class="flex items-center justify-between border-b py-2 px-4 bg-higher rounded-t-lg">
      <p class="pr-10 font-semibold">{{ title }}</p>
      <SLSearch
        class="bg-default max-w-[300px]"
        :placeholder="$t('experiment.index.config.table.search')"
        reverse
        v-model="search"
      />
    </div>
    <!-- table -->
    <div class="w-full">
      <!-- header -->
      <div class="line bg-higher border-b">
        <div class="header-item" v-for="item in column" :key="item.key" :title="item.title">{{ item.title }}</div>
      </div>
      <!-- body -->
      <div class="max-h-[500px] overflow-y-auto">
        <TableLine class="data-line" v-for="(line, index) in tableData" :key="line.key" :line="line" :index="index">
          <div class="body-item" v-for="item in line" :key="item.key" :title="item">
            {{ item }}
          </div>
        </TableLine>
      </div>
      <!-- empty body -->
      <div class="w-full py-4 flex justify-center" v-if="tableData.length === 0">
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
  @apply w-full rounded-lg outline;
  outline-color: var(--outline-default);
  outline-width: 1px;
}

.line {
  @apply w-full grid grid-cols-2 gap-3;
  grid-template-columns: 1fr 2fr;
}

.data-line {
  @apply line;
  &:last-child {
    @apply rounded-b-lg;
  }
}

.header-item {
  @apply border-r pl-4 my-3;
  &:last-child {
    @apply border-none pl-0;
  }
}

.body-item {
  @apply py-3 pl-4 pr-2 break-words;
  &:last-child {
    @apply border-none pl-0;
  }
}
</style>
