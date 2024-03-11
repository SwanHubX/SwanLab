<template>
  <!-- 这一层是为了保证表格宽度，设计到边框长度、背景色长度等 -->
  <div class="w-full" ref="wrapper">
    <!-- 表格 -->
    <div
      :class="{ 'gradient-last-line': lastRowGradient, 'table-border': tableBorder }"
      :style="{ width: typeof tableWidth === 'string' ? tableWidth : `${tableWidth}px` }"
      ref="table"
    >
      <!-- 表头 -->
      <div :class="{ 'sticky top-0 z-20': stickyHeader }">
        <div class="border-y flex bg-higher table-header">
          <!-- 表头项 -->
          <div
            v-for="(item, index) in column"
            :key="item.key"
            class="cell table-header-item"
            :class="activeColumnIndex === index ? 'bg-highest' : 'bg-higher'"
            @mouseover="handleMouseOver(index)"
            @mouseout="handleMouseOver(-1)"
            :style="{ width: elementWidths[index] }"
            ref="columnsRef"
          >
            <p
              class="overflow-hidden w-full h-full flex items-center"
              :class="item.style || 'px-2 py-3'"
              :title="item.title"
            >
              {{ item.title }}
            </p>
            <!-- 拖拽点 -->
            <button
              @mousedown="(e) => resize(e, index)"
              :class="{ 'button-hover-tip': hoverColumnIndex === index }"
              v-if="!column.unresizeable && !flexable"
            />
          </div>
        </div>
      </div>
      <!-- 表体, 按一行一行来渲染 -->
      <div
        v-for="(dataColumn, dataIndex) in data"
        :key="dataColumn"
        :class="{ 'hover:bg-higher': resizeIndex === -1 }"
        class="line"
      >
        <!-- 单元格 -->
        <div
          v-for="(item, index) in column"
          :key="item.key"
          :title="dataColumn[item.key]"
          class="cell flex items-center px-2 py-3"
          :class="[
            item.style,
            { 'hover:bg-primary-dimmest': resizeIndex === -1 },
            activeColumnIndex === index ? 'bg-higher' : 'bg-default'
          ]"
          @mouseover="handleMouseOver(index)"
          @mouseout="handleMouseOver(-1)"
          :style="{ width: elementWidths[index] }"
        >
          <div class="w-full overflow-hidden" v-if="item.slot">
            <slot :name="item.slot" v-bind:row="dataColumn" v-bind:index="dataIndex"></slot>
          </div>
          <!-- 文本格式 -->
          <div class="w-full overflow-hidden" v-else>{{ dataColumn[item.key] || '-' }}</div>
        </div>
      </div>
      <!-- 没有数据时，空占位 -->
      <div v-if="!data.length" class="flex justify-center py-3 border border-t-0">
        <div class="data-empty">{{ $t('common.table.empty') }}</div>
      </div>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 表格 —— 二次重构版
 * @file: ExprTable.vue
 * @since: 2023-12-28 22:37:21
 **/

import { computed } from 'vue'
import { onMounted } from 'vue'
import { ref, watch, inject } from 'vue'

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
  },
  // 最后一行渐变
  lastRowGradient: {
    type: Boolean,
    default: false
  },
  // 表格边框,除了上方一直存在以外，其他边框都需要设置tableBorder时才出现
  tableBorder: {
    type: Boolean,
    default: false
  },
  // 表头固定
  stickyHeader: {
    type: Boolean,
    default: false
  }
})

// ---------------------------------- 全局信息 ----------------------------------

const columnsRef = ref(null)
const widths = ref([])
const table = ref(null)
const wrapper = ref(null)

/**
 * 动态计算表格宽度
 * 需要动态计算表格宽度的原因是，overflow 出现滚动条之后，子元素的 100% 将不会计算超出部分
 * 这会导致一些样式问题，而动态计算之后，只需要将表格元素的父元素设置 overflow 即可
 */
const tableWidth = ref('100%')

// 宽度转化，直接控制每列宽度
const elementWidths = computed(() => {
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
// 当前被选中的列的索引
const activeColumnIndex = computed(() => {
  return resizeIndex.value === -1 ? hoverColumnIndex.value : resizeIndex.value
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

  // 监听视窗大小变化，以适配表格大小
  handleTableWith()
  window.addEventListener('resize', handleTableWith)
})

onUnmounted(() => {
  window.removeEventListener('resize', handleTableWith)
})

/**
 * 动态计算表格宽度
 */
const handleTableWith = () => {
  let width = 0
  widths.value.forEach(({ value }) => {
    width += value
  })
  if (wrapper.value?.offsetWidth > width) {
    if (wrapper.value.offsetWidth) tableWidth.value = wrapper.value.offsetWidth
  } else {
    tableWidth.value = width
  }
}

const layoutSiderBar = inject('isSideBarShow')
const siderWith = 288
watch(layoutSiderBar, (newV) => {
  if (props.flexable) return
  if (!newV) {
    tableWidth.value += siderWith
  } else {
    tableWidth.value -= siderWith
  }
})

// ---------------------------------- resize 相关 ----------------------------------

// 正在被重置宽度的列的索引
const resizeIndex = ref(-1)
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
  resizeIndex.value = index

  // 记录初始化数据
  startX.value = event.clientX
  if (elementWidths.value[index].endsWith('%')) {
    widths.value[index].value = columnsRef.value[index].offsetWidth
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
  widths.value[resizeIndex.value].value = newWidth
  columnsRef.value[resizeIndex.value] = newWidth
}

// 鼠标抬起，结束resize，清空状态，接触多余的事件绑定
const handleMouseup = () => {
  handleTableWith()
  document.removeEventListener('mousemove', handleMousemove)
  document.removeEventListener('mouseup', handleMouseup)
  document.body.style.userSelect = 'auto'
  document.body.style.cursor = ''
  startX.value = 0
  startWidth.value = 0
  resizeIndex.value = -1
}

// ---------------------------------- 处理鼠标移入移出事件，设置颜色 ----------------------------------
const handleMouseOver = (index) => {
  hoverColumnIndex.value = index
}
</script>

<style lang="scss" scoped>
.gradient-last-line {
  // 最后一行的cell添加伪元素渐变
  .line:last-child .cell:first-child {
    &:after {
      content: '';
      box-sizing: content-box;
      position: absolute;
      top: 100%;
      left: 0;
      height: 40px;
      width: 100%;
      border-right: 1px;
      border-style: solid;
      // 边框颜色从上到下渐变，从outline-default到transparent
      border-image: linear-gradient(to bottom, var(--outline-default), transparent) 1;
    }
  }
}

.table-border {
  .table-header {
    @apply border-x;
  }
  .line {
    @apply border-x;
  }
  .line:last-child {
    @apply border-b;
  }
}

.cell {
  @apply text-left whitespace-nowrap shrink-0;
  @apply box-border;
  height: 54px;

  &:first-child {
    @apply sticky left-0 z-10 border-r;
  }

  &:not(:first-child) {
    @apply z-0 relative;
  }
}

.line {
  @apply flex items-center;
  &:hover {
    .cell:not(:hover) {
      background-color: var(--background-higher);
    }
  }
}

.table-header {
  @apply relative;
  .table-header-item {
    &:hover {
      button {
        background-color: var(--foreground-dimmest);
      }
    }
    button {
      @apply w-1 h-full absolute right-0 top-0 cursor-col-resize;
      @apply hover:w-2 hover:bg-primary-dimmer active:bg-primary-dimmer;
      @apply transition-all duration-300;
    }
  }
}

.button-hover-tip {
  background-color: var(--foreground-dimmest);
}
</style>
