<template>
  <SLTable :header="header" :data="data" align="center" />
</template>

<script setup>
import SLTable from '@swanlab-vue/components/SLTable.vue'
import SLStatusLabel from '@swanlab-vue/components/SLStatusLabel.vue'
/**
 * @description: 表格组件使用示例，自己看着改
 * @file: TableExample.vue
 * @since: 2023-12-11 19:51:00
 **/
// 1x2的数组，代表表头
const header = ['Key', 'Value']
// Mx2的数组，代表表格内容
const data = [
  ['Key1', 'Value1'],
  ['Key2', { component: SLStatusLabel, props: { status: 1 } }],
  ['Key3', 'Value3']
]
</script>

<style lang="scss" scoped></style>
