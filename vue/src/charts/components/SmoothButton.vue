<template>
  <SLMenu class="right-0 w-96" corner-class="left-2.5" not-down>
    <div class="smooth-button">
      <SLIcon class="w-full h-full" icon="smooth" />
    </div>
    <template #pop>
      <div class="flex select-none flex-col w-full px-4 text-sm font-semibold gap-4">
        <div class="flex items-center gap-2">
          <span class="w-20 flex-shrink-0">{{ $t('common.chart.charts.share.smooth.smoothing') }}</span>
          <SLSlideBar v-model="_value" :max="1" :min="0" show-input @turn="handleTurn" @change="handleChange" />
        </div>
        <div class="flex items-center gap-2">
          <span class="w-20 flex-shrink-0">{{ $t('common.chart.charts.share.smooth.method') }}</span>
          <SLMenu menu-class="w-full" down down-position="top-2.5 right-2">
            <div class="w-full border p-2 text-left rounded">Gaussian</div>
            <template #pop>
              <div class="">Gaussian</div>
              <div class="">Mean</div>
              <div class="">Median</div>
            </template>
          </SLMenu>
        </div>
      </div>
    </template>
  </SLMenu>
</template>

<script setup>
/**
 * @description: 平滑按钮，包含菜单
 * @file: SmoothButton.vue
 * @since: 2024-03-04 21:13:50
 **/
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import SLSlideBar from '@swanlab-vue/components/SLSlideBar.vue'
import SLMenu from '@swanlab-vue/components/menu/SLMenu.vue'

const _value = ref(0)

// ---------------------------------- 处理slidebar的change和turn事件 ----------------------------------
const handleTurn = (direction, value) => {
  console.log('turn', direction, value)
}
const handleChange = (value) => {
  _value.value = Number(value)
}
</script>

<style lang="scss" scoped>
.smooth-button {
  @apply h-8 w-8 p-1 rounded text-dimmer hover:text-positive-default hover:bg-positive-dimmer;
}
</style>
