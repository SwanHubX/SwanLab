<template>
  <div class="w-full overflow-auto">
    <!-- 表格 -->
    <div class="border" ref="table" :style="{ width: table_width }">
      <!-- 表头 -->
      <div class="border-b flex bg-higher table-header">
        <!-- 表头项 -->
        <div
          v-for="(item, index) in column"
          :key="item.key"
          class="relative w-full cell shrink-0"
          :class="`${activeColumnIndex === index ? ' bg-highest' : 'bg-higher'} ${item.border ? 'border-r' : ''}`"
          @mouseover="() => (hoverColumnIndex = index)"
          @mouseout="() => (hoverColumnIndex = -1)"
          :style="{ width: element_widths[index] }"
          ref="columns"
        >
          <div
            class="overflow-hidden w-full h-full flex items-center"
            :class="item.style ? item.style : 'px-2 py-3'"
            :title="item.title"
          >
            {{ item.title }}
            <!-- 拖拽点 -->
            <span
              :class="`${activeColumnIndex === index ? 'bg-positive-dimmest' : ''} ${
                resize_index === index ? '!bg-primary-default' : ''
              }`"
              @mousedown="(e) => resize(e, index)"
              v-if="!column.unresizeable && !flexable"
            ></span>
          </div>
        </div>
      </div>
      <!-- 表体 -->
      <!-- 按一行一行来渲染 -->
      <div
        v-for="(dataColumn, dataIndex) in data"
        :key="dataColumn"
        :class="resize_index === -1 ? 'hover:bg-higher' : ''"
        class="flex items-center"
      >
        <!-- 单元格 -->
        <div
          v-for="(item, index) in column"
          :key="item.key"
          :title="dataColumn[item.key]"
          class="cell shrink-0 px-2 py-3"
          :class="`${'swanlab-table-column-' + index} ${resize_index === -1 ? 'hover:bg-primary-dimmest' : ''} ${
            item.border ? 'border-r ' : ''
          } ${item.style} ${activeColumnIndex === index ? 'bg-higher' : ''}`"
          @mouseover="() => (hoverColumnIndex = index)"
          @mouseout="() => (hoverColumnIndex = -1)"
          :style="{ width: element_widths[index] }"
        >
          <div v-if="item.slot">
            <slot :name="item.slot" v-bind:row="dataColumn" v-bind:index="dataIndex"></slot>
          </div>
          <!-- 文本格式 -->
          <div v-else>
            {{ dataColumn[item.key] || '-' }}
          </div>
        </div>
      </div>
      <!-- 没有数据时，空占位 -->
      <div v-if="!data.length" class="flex justify-center py-3">
        <div class="data-empty">{{ $t('common.table.empty') }}</div>
      </div>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 表格 —— 二次重构版
 * @file: SLTable.vue
 * @since: 2023-12-28 22:37:21
 **/

import { computed } from 'vue'
import { onMounted } from 'vue'
import { ref } from 'vue'

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
  },
  // 是否自动填充
  flexable: {
    type: Boolean,
    default: false
  }
})

// ---------------------------------- 全局信息 ----------------------------------

const columns = ref(null)
const widths = ref([])
const table = ref(null)

/**
 * 动态计算表格宽度
 * 需要动态计算表格宽度的原因是，overflow 出现滚动条之后，子元素的 100% 将不会计算超出部分
 * 这会导致一些样式问题，而动态计算之后，只需要将表格元素的父元素设置 overflow 即可
 */
const table_width = computed(() => {
  if (props.flexable) return '100%'
  let temp = 0
  widths.value.forEach(({ value }) => {
    temp += value
  })
  return temp + 2 + 'px'
})

// 宽度转化，直接控制每列宽度
const element_widths = computed(() => {
  if (props.flexable) {
    const len = props.column.length
    return new Array(len).fill((1 / len) * 100 + '%')
  }
  return widths.value.map((item) => {
    return item?.value ? item?.value + 'px' : ''
  })
})

// ---------------------------------- 样式相关 ----------------------------------

const hoverColumnIndex = ref(-1) // 被hover得列的索引
const activeColumnIndex = computed(() => {
  return resize_index.value === -1 ? hoverColumnIndex.value : resize_index.value
})

onMounted(() => {
  // 遍历列的设置
  props.column.forEach((column, index) => {
    const default_width = 150
    widths.value[index] = {
      value: 0,
      extended: false
    }
    // 获取所有的列宽
    // 如果自己设置了宽度
    if (column.width) {
      const target_width = column.width > default_width ? column.width : default_width
      widths.value[index].value = column.border ? target_width - 2 : target_width
      return
    }
    // 没自己设宽度，使用自然宽度
    widths.value[index].value = default_width
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
  if (element_widths.value[index].endsWith('%')) {
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
.cell {
  @apply overflow-hidden text-left whitespace-nowrap shrink-0;
}

.table-header span {
  @apply w-1.5 h-full absolute right-0 top-0 hover:bg-positive-higher hover:opacity-20 cursor-col-resize;
}
</style>
