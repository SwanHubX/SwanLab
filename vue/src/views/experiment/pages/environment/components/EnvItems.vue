<template>
  <div class="w-full border-b py-5 flex-col space-y-5" v-if="data?.length !== 0">
    <EnvItem
      :env-key="item.key"
      :env-value="item.value"
      v-for="item in data"
      :key="item.key"
      :highLight="item.highLight"
    />
  </div>
</template>

<script setup>
/**
 * @description: 环境分块 - 每个分块包含多条环境信息
 * @file: EnvItems.vue
 * @since: 2024-01-09 20:03:21
 **/

import EnvItem from './EnvItem.vue'

defineProps({
  data: {
    type: Array,
    default: () => {
      return []
    }
  }
})
</script>

<style lang="scss" scoped></style>
