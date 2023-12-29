<template>
  {{ tableWidth }}
  <div class="w-full relative">
    <div class="relative overflow-auto" :class="maxW">
      <div class="w-full">
        <table class="w-full border-collapse table-auto border overflow-auto">
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
                  <!-- 拖拽点 -->
                  <span
                    class="w-1.5 h-full absolute right-0 top-0 hover:bg-positive-dimmer hover:opacity-20 cursor-col-resize"
                    :class="hoverColumn === index ? 'bg-positive-highest' : ''"
                    @mousedown="(e) => resize(e, index)"
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
                :class="'swanlab-table-column-' + index"
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
  </div>
</template>

<script setup>
import { computed } from 'vue'
import { onMounted } from 'vue'
import { ref } from 'vue'

/**
 * @description: 表格 —— 二次重构版
 * @file: SLTable.vue
 * @since: 2023-12-28 22:37:21
 **/

const columns = ref(null)
const widths = ref([])
const tableWidth = computed(() => {
  let width = 0
  widths.value.forEach((item) => {
    width += item
  })
  return width
})

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
    // 获取所有的列宽
    widths.value[index] = columns.value[index].offsetWidth
    if (!column.width) return
    columns.value[index].setAttribute('width', column.width)
    widths.value[index] = column.width
  })
})

// ---------------------------------- resize 相关 ----------------------------------

const resize_index = ref(-1)
const startX = ref(0)
const startWidth = ref(0)
const resize = (event, index) => {
  console.log(`准备重置第${index}的宽度`)
  // 设置选取时不可以选中文本
  document.body.style.userSelect = 'none'
  document.body.style.cursor = 'col-resize'
  // 记录被改变大小的列的索引
  resize_index.value = index
  // 记录初始化数据
  console.log(columns.value[index])
  startX.value = event.clientX
  startWidth.value = columns.value[index].offsetWidth
  // 绑定全局事件
  document.addEventListener('mousemove', handleMousemove)
  document.addEventListener('mouseup', handleMouseup)
}

// 鼠标移动
const handleMousemove = (e) => {
  let newWidth = startWidth.value + (e.clientX - startX.value)
  widths.value[resize_index.value] = newWidth
  // 设置所有列
}

// 鼠标抬起，结束resize
const handleMouseup = () => {
  document.removeEventListener('mousemove', handleMousemove)
  document.removeEventListener('mouseup', handleMouseup)
  document.body.style.userSelect = 'auto'
  document.body.style.cursor = ''
}
</script>

<style lang="scss" scoped>
td,
th {
  @apply p-0 text-left whitespace-nowrap overflow-hidden;
}
</style>
