<template>
  <div class="gnip-table" ref="gnipTable">
    <div class="overflow-wrap w-full" :class="maxW" ref="overflowWrap" :style="{ width: initTableWidth + 'px' }">
      <div class="table-wrap w-full" :style="{ width: dynamicTableWidth + 'px' }">
        <table class="w-full">
          <!--  标签用于对表格中的列进行组合，以便对其进行格式化。 -->
          <colgroup>
            <col v-for="(item, index) in column" :key="index" ref="colgroupItems" />
          </colgroup>
          <!-- 表头 -->
          <thead>
            <tr>
              <th
                v-for="(item, index) in column"
                :key="index"
                class="gnip-th"
                :class="highlightColumnIndex === index ? ' bg-slate-200' : 'bg-[#f6f8fa]'"
                @mouseover="highlightColumn(highLight ? index : -1)"
                @mouseout="resetHighlight"
              >
                <span
                  class="block w-full whitespace-nowrap text-left"
                  :class="`${item.style || 'px-2'} ${item.fixed ? fixedTableWidth : ''}`"
                  :style="{
                    width: dragedRecord[index]
                      ? dragedRecord[index][0]
                        ? dragedRecord[index][1].toFixed(0) + 'px'
                        : ''
                      : columnDefault + 'px'
                  }"
                  >{{ item.title }}</span
                >
                <span
                  class="drag-line hover:!bg-positive-dimmer hover:opacity-20"
                  :class="`${activeDragIndex === index ? ' !bg-primary-dimmest' : ''} ${
                    highlightColumnIndex === index ? 'bg-positive-highest' : ''
                  } ${index + 1 === column.length ? 'pointer-events-none' : ''}`"
                  :style="{ background: index + 1 === column.length ? 'none' : '' }"
                  @mousedown="handleMouseDown(index, $event)"
                ></span>
              </th>
            </tr>
          </thead>
          <!-- 表体 -->
          <tbody v-if="data.length">
            <tr v-for="(dataColumn, dataIndex) in data" :key="dataColumn.id">
              <td
                v-for="(item, index) in column"
                :key="index"
                class="whitespace-nowrap overflow-hidden"
                :class="highlightColumnIndex === index ? 'hover:bg-primary-dimmest bg-[#ebf7ff]' : ''"
                @mouseover="highlightColumn(highLight ? index : -1)"
                @mouseout="resetHighlight"
              >
                <div
                  class="text-left overflow-hidden"
                  :class="`${item.style || 'px-2'}  ${item.fixed ? fixedTableWidth : ''}`"
                  :style="{
                    width: dragedRecord[index]
                      ? dragedRecord[index][0]
                        ? dragedRecord[index][1].toFixed(0) + 'px'
                        : ''
                      : columnDefault + 'px'
                  }"
                  v-if="item.slot"
                >
                  <slot :name="item.slot" v-bind:row="dataColumn" v-bind:index="dataIndex"></slot>
                </div>
                <div
                  class="text-left"
                  :class="`${item.style || 'px-2'}  ${item.fixed ? fixedTableWidth : ''}`"
                  :style="{
                    width: dragedRecord[index]
                      ? dragedRecord[index][0]
                        ? dragedRecord[index][1].toFixed(0) + 'px'
                        : ''
                      : columnDefault + 'px'
                  }"
                  v-else
                >
                  {{ dataColumn[item.key] || '-' }}
                </div>
              </td>
            </tr>
          </tbody>
          <tfoot v-if="!data.length">
            <tr>
              <td :colspan="column.length">
                <div class="data-empty">{{ $t('common.table.empty') }}</div>
              </td>
            </tr>
          </tfoot>
        </table>
        <!-- 拖拽线 -->
        <div class="drag-resize-line" :style="computedResizeLineStyle" v-if="showDragResizeLine"></div>
      </div>
    </div>
  </div>
</template>

<script>
/**
 * @description: 表格组件
 * 下面有组合式api的重构版，但是有一个小bug，即拖动单元格的时候，标准线一直在表格最左边
 * 并且，组合式版本有两种写法，目前先随便用一种，等重构完外部，再分析一下这里得优劣
 * @file: SLTable.vue
 * @since: 2023-12-24 15:30:00
 **/
export default {
  props: {
    // 表体展示数据
    data: {
      type: Array,
      default: () => {
        return []
      }
    },
    // 表头配置
    column: {
      type: Array,
      default: () => {
        return []
      }
    },
    // 最大宽度
    maxW: {
      type: String
    },
    // 是否高亮预览
    highLight: {
      type: Boolean,
      default: false
    }
  },
  data() {
    return {
      // 拖拽线左侧偏移量
      dragLineLeft: 0,
      // 当前拖拽的td索引项
      activeDragIndex: -1,
      // 开始拖拽那一刻，距离表格左侧的距离
      activeDragStartIndexOffsetLeft: 0,
      // 表格的宽度
      gnipTableScrollWidth: '',
      // 是否显示拖拽线
      showDragResizeLine: false,
      // 默认列宽为自动平分计算
      columnDefault: 0,
      // 记录对应列是被拖拽过
      dragedRecord: [],
      // 表格宽度，初始的时候固定宽度，超出滚动条，兼容多列数据放不下问题
      initTableWidth: 0,
      // 动态计算表格宽度
      dynamicTableWidth: 0,
      // 固定列样式
      fixedTableWidth: 'fixed  z-[2] bg-default',
      // 高亮行号
      highlightColumnIndex: -1
    }
  },
  computed: {
    computedResizeLineStyle() {
      return {
        left: this.dragLineLeft + 'px'
      }
    }
  },
  mounted() {
    this.init()
    window.addEventListener('resize', this.initTableScrollWidth)
  },
  methods: {
    // 初始化
    init() {
      this.initTableScrollWidth()
      // 初始化北拖拽的记录
      this.dragedRecord = new Array(this.column.length).fill(0).map((item, index) => {
        if (this.column[index].width) {
          // [是否被拖拽过或者设置过默认宽度，被拖拽的宽度或者具体拖拽后的宽度（实时更新)]
          return [true, this.column[index].width * 1]
        } else {
          return [false, 0]
        }
      })
    },
    // 初始化表格宽度
    initTableScrollWidth() {
      this.gnipTableScrollWidth = this.$refs.gnipTable.offsetWidth
      this.dynamicTableWidth = this.gnipTableScrollWidth - 1
      this.initTableWidth = this.gnipTableScrollWidth + 1
      // 自定义宽度的列数
      const autoWidthColumnCount = this.column.filter((item) => !item.width).length
      // 剩余可自定义宽度
      const resetAutoWidth = this.column.filter((item) => item.width).reduce((prev, now) => prev + now.width * 1, 0)
      // 剩余自适应宽度的列可分得的宽度，最低为0
      const columnResetWidth = (this.gnipTableScrollWidth - resetAutoWidth) / autoWidthColumnCount
      this.columnDefault = columnResetWidth > 0 ? columnResetWidth : 0
    },
    // 鼠标按下
    handleMouseDown(index, event) {
      // console.log("鼠标按下");
      // 计算拖拽线的距离
      this.computedDragLineOffsetStart(event.target)
      //记录当前拖拽的表头索引
      this.activeDragIndex = index
      // 显示拖拽线
      this.showDragResizeLine = true
      document.addEventListener('mousemove', this.handleMouseMove) //拖动中
      document.addEventListener('mouseup', this.handleMouseup) //抬起
    },
    // 拖动中
    handleMouseMove(event) {
      // console.log('拖动中')
      this.computedDragLineOffsetDraging(event)
      // 设置body鼠标样式为col-resize
      document.body.style.cursor = 'col-resize'
    },
    // 鼠标抬起
    handleMouseup(event) {
      this.computedDragLineOffsetDraging(event)
      // // 设置宽度
      this.setColGroupItemWidth()
      this.activeDragStartIndexOffsetLeft = 0
      // 隐藏拖拽线
      this.showDragResizeLine = false
      // 恢复body的鼠标样式
      document.body.style.cursor = ''

      // 移除事件
      document.removeEventListener('mousemove', this.handleMouseMove)
      document.removeEventListener('mousedown', this.handleMouseDown)
      document.removeEventListener('mouseup', this.handleMouseup)
    },
    // 计算拖拽线的偏移量
    computedDragLineOffsetStart(targetEle) {
      let parent = targetEle.offsetParent
      let left = targetEle.offsetLeft
      while (parent) {
        if (parent == this.$refs.overflowWrap) {
          // 说明找到了外层相对定位的父级
          break
        }
        left += parent.offsetLeft
        parent = parent.offsetParent
      }
      // 如果出现了水平滚动条，要加上水平滚动条卷走的宽度
      this.dragLineLeft = left + this.$refs.overflowWrap.scrollLeft
      // console.log(this.$refs.)
      // 记录按下刚开始拖拽的时候的位置
      !this.activeDragStartIndexOffsetLeft && (this.activeDragStartIndexOffsetLeft = left)
      return left
    },
    // 拖拽中计算偏移量
    computedDragLineOffsetDraging(event) {
      let clientX = event.clientX
      let targetEle = this.$refs.overflowWrap
      let x = targetEle.offsetLeft
      let parent = targetEle.offsetParent
      while (parent) {
        x += parent.offsetLeft
        parent = parent.offsetParent
      }
      let left = clientX - x < 0 ? 0 : clientX - x
      this.dragLineLeft = left + this.$refs.overflowWrap.scrollLeft
    },
    // 计算表头列的宽度
    setColGroupItemWidth() {
      // 需要增加的宽度
      let addWidth = this.dragLineLeft - this.activeDragStartIndexOffsetLeft

      this.$nextTick(() => {
        //  原来的宽度
        let { width } = this.$refs.colgroupItems[this.activeDragIndex].getBoundingClientRect()
        // 设置一个最小的宽度取值为30px，避免出现负值
        const computedWidth = addWidth + width < 70 ? 70 : addWidth + width
        // 动态设置表格宽度
        this.dynamicTableWidth =
          this.dynamicTableWidth + addWidth > this.initTableWidth
            ? (this.dynamicTableWidth += addWidth)
            : this.initTableWidth - 1
        // 标记记录当前列为拖拽过了
        this.dragedRecord[this.activeDragIndex] = [true, computedWidth]
        // 自定义宽度的列数(没被拖拽的重新计算列宽)
        // const autoWidthColumnCount = this.dragedRecord.filter((item) => !item[0]).length
        // 计算需要减去的宽度
        // const subtractWidth = this.dragedRecord
        //   .filter((item) => item[0])
        //   .reduce((prev, now) => {
        //     return now[1] + prev
        //   }, 0)
        for (let i = 0; i < this.$refs.colgroupItems.length; i++) {
          // 设置当前拖拽列的宽度
          if (i == this.activeDragIndex) {
            this.$refs.colgroupItems[this.activeDragIndex].setAttribute('width', computedWidth)
          } else {
            // 设置其它的列的宽度（排除配置项设置过的列)
            // const isComputed = (this.dragedRecord.find((item, index) => index == i) || [])[0]
            // !isComputed &&
            //   this.$refs.colgroupItems[i].setAttribute(
            //     'width',
            //     (this.dynamicTableWidth - subtractWidth) / autoWidthColumnCount
            //   )
          }
        }
        // 重置表头索引
        this.activeDragIndex = -1
        // 重置拖拽数据
        // this.dragLineLeft = 0
      })
    },
    // 设置行高亮
    highlightColumn(index) {
      this.highlightColumnIndex = index
    },
    // 删除高亮
    resetHighlight() {
      this.highlightColumnIndex = -1
    }
  }
}
</script>

<!-- <script setup>
import { ref, onMounted, nextTick } from 'vue'

const props = defineProps({
  data: {
    type: Array,
    default: () => {
      return []
    }
  },
  column: {
    type: Array,
    default: () => {
      return []
    }
  }
})

// 左拖拽偏移量
const dragLineLeft = ref(0)
// 当前拖拽的td索引项
const activeDragIndex = ref(-1)
// 开始拖拽那一刻，距离表格左侧
const activeDragStartIndexOffsetLeft = ref(0)
// 表格的宽度
const gnipTableScrollWidth = ref('')
// 是否显示拖线
const showDragResizeLine = ref(false)
// 默认列宽为自动平分计算
const columnDefault = ref(0)
// 记录对应列是被拖拽过
const dragedRecord = ref([])
// 表格宽度，初始的时候固定宽度，超出滚动条，兼容多列数据放不下问题
const initTableWidth = ref(0)
// 动态计算表格宽度
const dynamicTableWidth = ref(0)

// 计算拖拽后的样式
const computedResizeLineStyle = ref({
  left: dragLineLeft.value + 'px'
})

// 初始化
const init = () => {
  initTableScrollWidth()
  // 初始化北拖拽的记录
  dragedRecord.value = new Array(props.column.length).fill(0).map((item, index) => {
    if (props.column[index].width) {
      // [是否被拖拽过或者设置过默认宽度，被拖拽的宽度或者具体拖拽后的宽度（实时更新)]
      return [true, props.column[index].width * 1]
    } else {
      return [false, 0]
    }
  })
}

const gnipTable = ref(null)
// 初始化表格宽度
const initTableScrollWidth = () => {
  gnipTableScrollWidth.value = gnipTable.value.offsetWidth
  dynamicTableWidth.value = gnipTableScrollWidth.value - 1
  initTableWidth.value = gnipTableScrollWidth.value + 1
  // 自定义宽度的列数
  const autoWidthColumnCount = props.column.filter((item) => !item.width).length
  // 剩余可自定义宽度
  const resetAutoWidth = props.column.filter((item) => item.width).reduce((prev, now) => prev + now.width * 1, 0)
  // 剩余自适应宽度的列可分得的宽度，最低为0
  const columnResetWidth = (gnipTableScrollWidth.value - resetAutoWidth) / autoWidthColumnCount
  columnDefault.value = columnResetWidth > 0 ? columnResetWidth : 0
}

// 鼠标按下
const handleMouseDown = (index, event) => {
  // 计算拖拽线的距离
  computedDragLineOffsetStart(event.target)
  //记录当前拖拽的表头索引
  activeDragIndex.value = index
  // 显示拖拽线
  showDragResizeLine.value = true
  document.addEventListener('mousemove', handleMouseMove) //拖动中
  document.addEventListener('mouseup', handleMouseup) //抬起
}

// 拖动中
const handleMouseMove = (event) => {
  // console.log("拖动中");
  computedDragLineOffsetDraging(event)
  // 设置body鼠标样式为col-resize
  document.body.style.cursor = 'col-resize'
}

// 鼠标抬起
const handleMouseup = (event) => {
  computedDragLineOffsetDraging(event)
  // 设置宽度
  setColGroupItemWidth()
  activeDragStartIndexOffsetLeft.value = 0
  // 隐藏拖拽线
  showDragResizeLine.value = false
  // 恢复body的鼠标样式
  document.body.style.cursor = ''

  // 移除事件
  document.removeEventListener('mousemove', handleMouseMove)
  document.removeEventListener('mousedown', handleMouseDown)
  document.removeEventListener('mouseup', handleMouseup)
}

const overflowWrap = ref(null)
// 计算拖拽线的偏移量
const computedDragLineOffsetStart = (targetEle) => {
  let parent = targetEle.offsetParent
  let left = targetEle.offsetLeft

  while (parent) {
    if (parent == overflowWrap.value) {
      // 说明找到了外层相对定位的父级
      break
    }
    left += parent.offsetLeft
    parent = parent.offsetParent
  }
  // 如果出现了水平滚动条，要加上水平滚动条卷走的宽度
  dragLineLeft.value = left + overflowWrap.value.scrollLeft
  // 记录按下刚开始拖拽的时候的位置
  !activeDragStartIndexOffsetLeft.value && (activeDragStartIndexOffsetLeft.value = left)

  return left
}

// 拖拽中计算偏移量
const computedDragLineOffsetDraging = (event) => {
  let clientX = event.clientX
  let targetEle = overflowWrap.value
  let x = targetEle.offsetLeft
  let parent = targetEle.offsetParent

  while (parent) {
    x += parent.offsetLeft
    parent = parent.offsetParent
  }

  let left = clientX - x < 0 ? 0 : clientX - x
  dragLineLeft.value = left + overflowWrap.value.scrollLeft
}

const colgroupItems = ref(null)
// 计算表头列的宽度
const setColGroupItemWidth = () => {
  // 需要增加的宽度
  let addWidth = dragLineLeft.value - activeDragStartIndexOffsetLeft.value

  nextTick(() => {
    //  原来的宽度
    let { width } = colgroupItems.value[activeDragIndex.value].getBoundingClientRect()
    // 设置一个最小的宽度取值为30px，避免出现负值
    const computedWidth = addWidth + width < 30 ? 30 : addWidth + width

    // 动态设置表格宽度
    dynamicTableWidth.value =
      dynamicTableWidth.value + addWidth > initTableWidth.value
        ? (dynamicTableWidth.value += addWidth)
        : initTableWidth.value - 1

    // 标记记录当前列为拖拽过了
    dragedRecord.value[activeDragIndex.value] = [true, computedWidth]
    // 自定义宽度的列数(没被拖拽的重新计算列宽)
    const autoWidthColumnCount = dragedRecord.value.filter((item) => !item[0]).length
    // 计算需要减去的宽度
    const subtractWidth = dragedRecord.value.filter((item) => item[0]).reduce((prev, now) => now[1] + prev, 0)

    for (let i = 0; i < colgroupItems.value.length; i++) {
      // 设置当前拖拽列的宽度
      if (i == activeDragIndex.value) {
        colgroupItems.value[i].setAttribute('width', computedWidth)
      } else {
        // 设置其它的列的宽度（排除配置项设置过的列)
        const isComputed = (dragedRecord.value.find((item, index) => index == i) || [])[0]
        !isComputed &&
          colgroupItems.value[i].setAttribute('width', (dynamicTableWidth.value - subtractWidth) / autoWidthColumnCount)
      }
    }
  })
}

onMounted(() => {
  init()
})
</script> -->

<style lang="scss">
.gnip-table {
  position: relative;
  // .table-wrap {
  //   overflow: auto;
  // }
  .overflow-wrap {
    position: relative;
    overflow: auto;
  }
  table {
    border-collapse: collapse;
    table-layout: auto;
    border: 1px solid #e8eaec;
    .data-empty {
      text-align: center;
    }
    .gnip-th {
      position: relative;
      // background-color: #f6f8fa;
      padding: 8px 0;
      &:hover span:last-child {
        @apply bg-positive-highest;
      }
      .drag-line {
        position: absolute;
        width: 5px;
        height: 100%;
        right: 0;
        top: 0;
        cursor: col-resize;
        user-select: none;
        z-index: 1;
      }
    }
  }
  table,
  th,
  td {
    border: 1px solid #e8eaec;
    text-align: center;
    word-break: break-all;
  }
  thead {
    .th {
      background-color: #f8f8f9;
    }
  }
  .drag-resize-line {
    height: 100%;
    width: 1px;
    border-right: 1px dashed green;
    position: absolute;
    left: 0;
    top: 0;
    z-index: 100;
  }
  tbody {
    tr {
      &:hover {
        background: #ebf7ff;
      }
    }
    td {
      height: 48px;
      box-sizing: border-box;
    }
  }
}
</style>
