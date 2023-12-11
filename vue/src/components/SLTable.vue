<template>
  <div class="border rounded overflow-hidden">
    <table class="sl-table" :class="[alignClass]">
      <tr>
        <th v-for="item in header" :key="item">
          <component :is="item.component" :props="item.props" v-if="shouldBeComponent(item)" />
          <template v-else>
            {{ item }}
          </template>
        </th>
      </tr>
      <!-- 数据部分 -->
      <tr v-for="d in data" :key="d[0]">
        <td v-for="item in d" :key="item">
          <component :is="item.component" v-bind="item.props" v-if="shouldBeComponent(item)" />
          <template v-else>
            {{ item }}
          </template>
        </td>
      </tr>
    </table>
  </div>
</template>

<script setup>
/**
 * @description: 普通的表格组件，输入data数据，输出表格;目前表格每列宽度固定
 * @param { Array } header 表头数据，格式为1xN的数组，N代表表头的数量
 * @param { Array } data 表格数据，格式为MxN的数组，M代表表格的行数，N代表表格的列数,N与header的长度相同
 *
 * 无论是data还是header，里面的元素支持以组件的形式传入，如果不是组件，会自动转换为字符串填充
 * 如果是组件，改为对象形式，格式为{component: 组件对象, props: 组件的props}
 *
 * @param { String } align 表格的对齐方式，可选值为left、center、right，默认为left
 * @file: SLTable.vue
 * @since: 2023-12-11 19:45:01
 **/

import { computed } from 'vue'
const props = defineProps({
  header: {
    type: Array,
    required: true
  },
  data: {
    type: Array,
    required: true
  },
  align: {
    type: String,
    default: 'left'
  }
})

const alignClass = computed(() => {
  return 'text-' + props.align
})

// ---------------------------------- 动态组件与字符串切换----------------------------------
/**
 * 切换为组件渲染的标志是这是一个对象且其中包含component字段
 */
const shouldBeComponent = (o) => {
  return typeof o === 'object' && !!o.component
}
</script>

<style lang="scss" scoped>
.sl-table {
  @apply w-full;
  tr {
    &:first-child {
      @apply bg-dimmer;
    }
    &:not(:first-child) {
      @apply border-t;
    }
  }

  th {
    text-align: inherit;
  }

  th,
  td {
    @apply px-5 py-2.5;

    &:not(:last-child) {
      @apply border-r;
    }
  }
}
</style>
