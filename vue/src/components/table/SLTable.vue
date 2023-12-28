<template>
  <div class="w-full relative">
    <div class="relative overflow-auto" :class="maxW">
      <table class="w-full border-collapse table-auto border">
        <!-- 标签用于对表格中的列进行组合，以便对其进行格式化 -->
        <colgroup>
          <col v-for="(item, index) in column" :key="item.key + item.slot" />
        </colgroup>
        <!-- 表头 -->
        <thead>
          <tr>
            <th v-for="(item, index) in column" :key="index" class="py-2">
              <span>{{ item.title }}</span>
              <!-- 拖拽点 -->
              <span></span>
            </th>
          </tr>
        </thead>
        <tbody></tbody>
      </table>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 表格 —— 二次重构版
 * @file: SLTable.vue
 * @since: 2023-12-28 22:37:21
 **/

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
</script>

<style lang="scss" scoped>
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
