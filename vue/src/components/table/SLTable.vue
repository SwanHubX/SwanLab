<template>
  <div class="w-full relative">
    <div class="relative overflow-auto" :class="maxW">
      <table class="w-full border-collapse table-auto border">
        <!-- 标签用于对表格中的列进行组合，以便对其进行格式化 -->
        <colgroup>
          <col
            v-for="(item, index) in column"
            :key="item.key + item.slot + index"
            :class="`${hoverColumn === index ? activeColumnBackground : ''} ${'swanlab-table-column-' + index}`"
            ref="columns"
          />
        </colgroup>
        <!-- 表头 -->
        <thead class="border-bottom">
          <tr>
            <th
              v-for="(item, index) in column"
              :key="item.key"
              class="relative overflow-hidden"
              :class="hoverColumn === index ? 'bg-slate-200' : 'bg-slate-50'"
              @mouseover="() => (hoverColumn = index)"
              @mouseout="() => (hoverColumn = -1)"
            >
              <div class="overflow-hidden" :class="item.style ? item.style : 'px-2 py-3'">
                {{ item.title }}
                <span
                  class="w-1.5 h-full absolute right-0 top-0 hover:bg-positive-dimmer hover:opacity-20 cursor-col-resize"
                  :class="hoverColumn === index ? 'bg-positive-highest' : ''"
                ></span>
              </div>
            </th>
          </tr>
        </thead>
        <!-- 表体 -->
        <tbody>
          <!-- 每一行 -->
          <tr v-for="(dataColumn, dataIndex) in data" :key="dataColumn" class="hover:bg-blue-50">
            <!-- 单元格 -->
            <td
              v-for="(item, index) in column"
              :key="item.key"
              class="hover:bg-blue-100 overflow-hidden"
              @mouseover="() => (hoverColumn = index)"
              @mouseout="() => (hoverColumn = -1)"
            >
              <div class="overflow-hidden" :class="item.style ? item.style : 'px-2 py-3'">
                <div v-if="item.slot">
                  <slot :name="item.slot" v-bind:row="dataColumn" v-bind:index="dataIndex"></slot>
                </div>
                <!-- 文本格式 -->
                <div v-else>
                  {{ dataColumn[item.key] || '-' }}
                </div>
              </div>
            </td>
          </tr>
        </tbody>
      </table>
    </div>
  </div>
</template>

<script setup>
import { onMounted } from 'vue'
import { ref } from 'vue'

/**
 * @description: 表格 —— 二次重构版
 * @file: SLTable.vue
 * @since: 2023-12-28 22:37:21
 **/

const columns = ref(null)

// ---------------------------------- 组件接口 ----------------------------------

const props = defineProps({
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
  },
  // 表格最大宽度
  maxW: {
    type: String,
    default: 'none'
  },
  // 是否高亮预览
  highLight: {
    type: Boolean,
    default: true
  }
})

// ---------------------------------- 样式相关 ----------------------------------

const activeColumnBackground = 'bg-blue-50'
const hoverColumn = ref(-1) // 被hover得列的索引
onMounted(() => {
  // 遍历列的设置
  props.column.forEach((column, index) => {
    if (!column.width) return
    columns.value[index].setAttribute('width', column.width)
  })
})

// ---------------------------------- resize 相关 ----------------------------------
</script>

<style lang="scss" scoped>
td,
th {
  @apply p-0 text-left whitespace-nowrap overflow-hidden;
}
</style>
