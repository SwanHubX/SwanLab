<template>
  <div class="gnip-table" ref="gnipTable">
    <div class="overflow-wrap" ref="overflowWrap" :style="{ width: initTableWidth + 'px' }">
      <div class="table-wrap" :style="{ width: dynamicTableWidth + 'px' }">
        <table :style="{ width: dynamicTableWidth + 'px' }">
          <!--  标签用于对表格中的列进行组合，以便对其进行格式化。 -->
          <colgroup>
            <col
              :width="item.width || columnDefault"
              v-for="(item, index) in column"
              :key="index"
              ref="colgroupItems"
            />
          </colgroup>
          <!-- 表头 -->
          <thead>
            <tr>
              <th v-for="(item, index) in column" :key="index" class="gnip-th">
                <span class="cell-title">{{ item.title }}</span>
                <span class="drag-line" @mousedown="handleMouseDown(index, $event)"></span>
              </th>
            </tr>
          </thead>
          <!-- 表体 -->
          <tbody v-if="data.length">
            <tr v-for="(dataColumn, dataIndex) in data" :key="dataColumn.id">
              <td v-for="(item, index) in column" :key="index">
                <div class="content-cel" v-if="item.slot">
                  <slot :name="item.slot" v-bind:row="dataColumn" v-bind:index="dataIndex"></slot>
                </div>
                <div class="content-cel" v-else>
                  {{ dataColumn[item.key] }}
                </div>
              </td>
            </tr>
          </tbody>
          <tfoot v-if="!data.length">
            <tr>
              <td :colspan="column.length">
                <div class="data-empty">暂无数据</div>
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
export default {
  props: {
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
      dynamicTableWidth: 0
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
      // console.log("拖动中");
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
        const computedWidth = addWidth + width < 30 ? 30 : addWidth + width
        // 动态设置表格宽度
        this.dynamicTableWidth =
          this.dynamicTableWidth + addWidth > this.initTableWidth
            ? (this.dynamicTableWidth += addWidth)
            : this.initTableWidth - 1
        // 标记记录当前列为拖拽过了
        this.dragedRecord[this.activeDragIndex] = [true, computedWidth]
        // 自定义宽度的列数(没被拖拽的重新计算列宽)
        const autoWidthColumnCount = this.dragedRecord.filter((item) => !item[0]).length
        // 计算需要减去的宽度
        const subtractWidth = this.dragedRecord
          .filter((item) => item[0])
          .reduce((prev, now) => {
            return now[1] + prev
          }, 0)
        for (let i = 0; i < this.$refs.colgroupItems.length; i++) {
          // 设置当前拖拽列的宽度
          if (i == this.activeDragIndex) {
            this.$refs.colgroupItems[this.activeDragIndex].setAttribute('width', computedWidth)
          } else {
            // 设置其它的列的宽度（排除配置项设置过的列)
            const isComputed = (this.dragedRecord.find((item, index) => index == i) || [])[0]
            !isComputed &&
              this.$refs.colgroupItems[i].setAttribute(
                'width',
                (this.dynamicTableWidth - subtractWidth) / autoWidthColumnCount
              )
          }
        }
      })
    }
  }
}
</script>

<!-- <script setup>
/**
 * @description: 表格组件
 * @file: SLTable.vue
 * @since: 2023-12-24 15:30:00
 **/
import { ref } from 'vue'

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
</script> -->

<style lang="scss">
.gnip-table {
  position: relative;
  .table-wrap {
    // overflow: auto;
  }
  .overflow-wrap {
    position: relative;
    overflow: auto;
  }
  table {
    border-collapse: collapse;
    table-layout: fixed;
    border: 1px solid #e8eaec;
    .data-empty {
      text-align: center;
    }
    .gnip-th {
      position: relative;
      background-color: #f8f8f9;
      padding: 8px 0;
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
