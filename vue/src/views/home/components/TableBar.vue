<template>
  <div class="w-full flex justify-between items-center">
    <SLSearch class="max-w-96" :placeholder="$t('experiment.index.header.table-bar.placeholder')"></SLSearch>
    <div class="flex gap-5 items-center pl-10">
      <SLCheck :label="$t('experiment.index.header.table-bar.summary')" />
      <SLButton
        theme="default"
        class="px-3 py-2 rounded-lg flex items-center gap-2"
        @click="downloadCsv(tableHead, tableBody)"
      >
        <SLIcon icon="home" class="w-4 h-4" />
        <span>{{ $t('experiment.index.header.table-bar.export') }}</span>
      </SLButton>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 表格上方功能栏：搜索、筛选、导出
 * @file: TableBar.vue
 * @since: 2024-01-16 17:52:07
 **/
import SLCheck from '@swanlab-vue/components/SLCheck.vue'
import { formatTime } from '@swanlab-vue/utils/time'
import { t } from '@swanlab-vue/i18n'

defineProps({
  tableBody: {
    type: Array,
    default: () => []
  },
  tableHead: {
    type: Array,
    default: () => []
  }
})

/**
 * @description 纯前端实现将表格数据导出为csv格式文件
 * @param {Array} headers 表格头配置项，headers中的key值与data中每一个item的属性名一一对应
 * @param {Array} data 表格数据
 * @param {String} fileName 导出的文件名称
 */
function downloadCsv(header, data, fileName = '导出结果.csv') {
  console.log(data)
  if (!header || !data || !Array.isArray(header) || !Array.isArray(data) || !header.length || !data.length) {
    return
  }
  var csvContent = 'data:text/csv;charset=utf-8,\ufeff'
  // 头部数据，以逗号分隔
  const _header = header.map((h) => h.title).join(',')
  csvContent += _header + '\n'

  // 遍历表体的每一行
  data.forEach((item, index) => {
    let dataString = ''
    // 遍历一行中的每一列
    for (let i = 0; i < header.length; i++) {
      const head = header[i]
      const label = head.key || head.slot
      dataString += (head.type ? checkType(item, head.type, label) : item[label] || '-') + ','
    }
    csvContent += index < data.length ? dataString.replace(/,$/, '\n') : dataString.replace(/,$/, '')
  })

  const a = document.createElement('a')
  a.href = encodeURI(csvContent)
  a.download = fileName
  a.click()
  window.URL.revokeObjectURL(csvContent)
}

/**
 * 特殊情况，即有些数据是需要处理，或具有特殊层级
 * @param {obj} item 对应行的数据对象
 * @param {string} type 类型
 * @param {string} key key 或 slot
 */
const checkType = (item, type, key) => {
  switch (type) {
    case 'config':
      return item.config[key]
    case 'create_time':
      return formatTime(item.create_time)
    case 'status':
      return t(`experiment.status.${item.status}`)
  }
}
</script>

<style lang="scss" scoped></style>