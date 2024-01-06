<template>
  <div class="w-full relative">
    <div class="relative w-full" :class="maxW">
      <div class="w-full overflow-auto border" :class="data.length ? '' : 'overflow-hidden'">
        <table class="w-full">
          <!-- 标签用于对表格中的列进行组合，以便对其进行格式化 -->
          <colgroup>
            <col
              v-for="(item, index) in column"
              :key="item.key + item.slot + index"
              :class="`${activeColumnIndex === index ? 'bg-higher' : ''} ${'swanlab-table-column-' + index}`"
              ref="columns"
            />
          </colgroup>
          <!-- 表头 -->
          <thead class="border-bottom">
            <tr class="border-b">
              <th
                v-for="(item, index) in column"
                :key="item.key"
                class="relative overflow-hidden"
                :class="`${activeColumnIndex === index ? ' bg-highest' : 'bg-higher'} ${item.border ? 'border-r' : ''}`"
                @mouseover="() => (hoverColumnIndex = index)"
                @mouseout="() => (hoverColumnIndex = -1)"
              >
                <div
                  class="overflow-hidden"
                  :class="item.style ? item.style : 'px-2 py-3'"
                  :title="item.title"
                  :style="{ width: widths[index]?.value + 'px' }"
                >
                  {{ item.title }}
                  <!-- 拖拽点 -->
                  <span
                    class="w-1.5 h-full absolute right-0 top-0 hover:bg-positive-higher hover:opacity-20 cursor-col-resize"
                    :class="`${activeColumnIndex === index ? 'bg-positive-dimmest' : ''} ${
                      resize_index === index ? '!bg-primary-default' : ''
                    }`"
                    @mousedown="(e) => resize(e, index)"
                    v-if="!column.unresizeable"
                  ></span>
                </div>
              </th>
            </tr>
          </thead>
          <!-- 表体 -->
          <tbody>
            <!-- 每一行 -->
            <tr
              v-for="(dataColumn, dataIndex) in data"
              :key="dataColumn"
              :class="resize_index === -1 ? 'hover:bg-higher' : ''"
            >
              <!-- 单元格 -->
              <td
                v-for="(item, index) in column"
                :key="item.key"
                :title="dataColumn[item.key]"
                class="overflow-hidden"
                :class="`${'swanlab-table-column-' + index} ${resize_index === -1 ? 'hover:bg-primary-dimmest' : ''} ${
                  item.border ? 'border-r' : ''
                }`"
                @mouseover="() => (hoverColumnIndex = index)"
                @mouseout="() => (hoverColumnIndex = -1)"
              >
                <div :class="item.style ? item.style : 'px-2 py-3'" :style="{ width: widths[index]?.value + 'px' }">
                  <div v-if="item.slot" :style="{ width: widths[index]?.value + 'px' }">
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
        <!-- 没有数据时，空占位 -->
        <div v-if="!data.length" class="flex justify-center py-3">
          <div class="data-empty">{{ $t('common.table.empty') }}</div>
        </div>
      </div>
    </div>
  </div>
  {{ widths }}
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

const hoverColumnIndex = ref(-1) // 被hover得列的索引
const activeColumnIndex = computed(() => {
  return resize_index.value === -1 ? hoverColumnIndex.value : resize_index.value
})

onMounted(() => {
  // 遍历列的设置
  props.column.forEach((column, index) => {
    const auto_width = columns.value[index].offsetWidth
    widths.value[index] = {
      value: 0,
      extended: false
    }
    // 获取所有的列宽
    // 如果自己设置了宽度
    if (column.width) {
      const target_width = column.width > auto_width ? column.width : auto_width
      widths.value[index].value = column.border ? target_width - 2 : target_width
      return
    }
    // 没自己设宽度，使用自然宽度
    widths.value[index].value = auto_width
  })
})

// ---------------------------------- resize 相关 ----------------------------------

// 正在被重置宽度的列的索引
const resize_index = ref(-1)
// 开始时鼠标位置
const startX = ref(0)
// 开始时，目标列的宽度
const startWidth = ref(0)
// 最小的列宽
const minCellWidth = ref(70)

/**
 * 准备重置宽度
 * @param {object} event 点击事件对象
 * @param {number} index 目标列的索引
 */
const resize = (event, index) => {
  console.log(`准备重置第${index}的宽度`)
  // 设置选取时不可以选中文本
  document.body.style.userSelect = 'none'
  document.body.style.cursor = 'col-resize'

  // 记录被改变大小的列的索引
  resize_index.value = index

  // 记录初始化数据
  startX.value = event.clientX
  if (!widths.value[index].value) {
    widths.value[index].value = columns.value[index].offsetWidth
  }
  startWidth.value = widths.value[index].value

  // 绑定全局事件
  document.addEventListener('mousemove', handleMousemove)
  document.addEventListener('mouseup', handleMouseup)
}

// 鼠标移动
const handleMousemove = (e) => {
  // 计算新宽度
  let newWidth = startWidth.value + (e.clientX - startX.value)
  // 对新宽度设置最小值
  newWidth = newWidth < minCellWidth.value ? minCellWidth.value : newWidth
  // 重置宽度
  widths.value[resize_index.value].value = newWidth
  columns.value[resize_index.value] = newWidth
}

// 鼠标抬起，结束resize，清空状态，接触多余的事件绑定
const handleMouseup = () => {
  document.removeEventListener('mousemove', handleMousemove)
  document.removeEventListener('mouseup', handleMouseup)
  document.body.style.userSelect = 'auto'
  document.body.style.cursor = ''
  startX.value = 0
  startWidth.value = 0
  resize_index.value = -1
}
</script>

<style lang="scss" scoped>
td,
th {
  @apply p-0 text-left whitespace-nowrap overflow-hidden;
}
</style>
