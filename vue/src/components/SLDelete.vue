<template>
  <SLButton
    @click="click"
    theme="negative"
    hollow
    class="px-3 py-1.5 rounded-lg"
    :disabled-tip="$t('common.delete.not-allowed')"
  >
    <div class="text-sm font-semibold flex items-center gap-2">
      <SLIcon icon="trash" class="w-4 h-4" />
      {{ $t('common.delete.button') }}
    </div>
  </SLButton>
</template>

<script setup>
/**
 * @description: 删除套件
 * @file: SLDelete.vue
 * @since: 2024-01-04 14:27:27
 *
 * 套件主要含有三部分：
 * 1. 适配于 type 的描述部分：不同的删除弹窗大概率会有不同的标题和内容
 * 2. 确认删除弹窗
 * 3. 取消/确认按钮的回调
 **/

import SLButton from './SLButton.vue'
import SLIcon from './SLIcon.vue'
import { confirm } from './comfirm'
import { t } from '@swanlab-vue/i18n'

const props = defineProps({
  type: String,
  validator: (value) => {
    // 目前type只适配了project和experiment两种模式
    return ['project', 'experiment'].includes(value)
  }
})

const emits = defineEmits(['confirm', 'cancel'])

const click = () => {
  confirm(t(`common.delete.title.${props.type}`), t(`common.delete.content.${props.type}`), {}).then(() => {
    emits('confirm')
  })
}
</script>

<style lang="scss" scoped></style>
