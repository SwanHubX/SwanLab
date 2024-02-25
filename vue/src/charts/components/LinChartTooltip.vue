<template>
  <div class="lc-tooltip" ref="toolTipRef" v-show="isShow" :style="{ width: getTooltipWidth() }" :key="key">
    <div class="lc-tooltip-item-zoom" v-if="detail">
      <p class="lc-tooltip-color"></p>
      <p class="lc-tooltip-step">{{ $t('common.chart.charts.share.step') }}</p>
      <p class="lc-tooltip-value">{{ $t('common.chart.charts.share.value') }}</p>
      <p class="lc-tooltip-time">{{ $t('common.chart.charts.share.time') }}</p>
      <p class="lc-tooltip-tag">{{ $t('common.chart.charts.share.tag') }}</p>
    </div>
    <template v-if="detail && items.length">
      <div class="lc-tooltip-item-zoom" v-for="item in items" :key="item.color" :style="{ color: item.color }">
        <!-- 颜色 -->
        <span class="lc-tooltip-color lc-tooltip-color-rect"></span>
        <!-- 步数 -->
        <span class="lc-tooltip-step">{{ item.data.index }}</span>
        <!-- 数据 -->
        <span class="lc-tooltip-value">{{ formatNumber2SN(item.data.data) }}</span>
        <!-- 时间 -->
        <span class="lc-tooltip-time">{{ formatTime(item.data.create_time) }}</span>
        <!-- 标签 -->
        <span class="lc-tooltip-tag">{{ item.data.series }}</span>
      </div>
    </template>
    <template v-else-if="items.length">
      <div class="lc-tooltip-item-no-zoom" v-for="item in items" :key="item.data.color" :style="{ color: item.color }">
        <!-- 颜色 -->
        <span class="lc-tooltip-color lc-tooltip-color-rect"></span>
        <!-- 步数 -->
        <span class="lc-tooltip-step">{{ item.data.index }}</span>
        <!-- 数据 -->
        <span class="lc-tooltip-value">{{ formatNumber2SN(item.data.data) }}</span>
        <!-- 标签 -->
        <span class="lc-tooltip-tag">{{ item.data.series }}</span>
      </div>
    </template>
    <p class="lc-tooltip-tip">{{ tip }}</p>
  </div>
</template>

<script setup>
/**
 * @description: 自定义折线图提示框
 * @file: LinChartTooltip.vue
 * @since: 2024-02-24 15:36:59
 **/
import { reactive, ref, inject } from 'vue'
import { isApple } from '@swanlab-vue/utils/browser'
import { t } from '@swanlab-vue/i18n'
const props = defineProps({
  // 是否显示详细版
  detail: {
    type: Boolean,
    default: false
  }
})

const tip = isApple ? t('common.chart.charts.line.copy.apple') : t('common.chart.charts.line.copy.windows')
const tooltipWidth = props.detail ? 400 : 256

const getTooltipWidth = () => {
  return `${tooltipWidth}px`
}
// 提示框数据
const items = ref([])
// 显示模式，分为详细版和简单版
const toolTipRef = ref(null)
const isShow = ref(false)
const tooltipXOffset = 50
const key = ref(null)
/**
 * 显示提示框，传入父元素宽度和显示x轴位置，计算提示框位置
 */
const show = (data, width, x) => {
  // console.log('show', data, width, x)
  // key.value = Math.random()
  items.value = data
  isShow.value = true
  const left = parseFloat(x.split('px')[0])
  // console.log(left + tooltipWidth, width)
  if (left + tooltipWidth > width) {
    toolTipRef.value.style.right = `${width - left + tooltipXOffset}px`
    toolTipRef.value.style.left = 'auto'
  } else {
    toolTipRef.value.style.left = `${left + tooltipXOffset}px`
    toolTipRef.value.style.right = 'auto'
  }
}

/**
 * 隐藏提示框，直接隐藏
 */
const hide = () => {
  isShow.value = false
}

// ---------------------------------- 依赖获取 ----------------------------------
const formatNumber2SN = inject('formatNumber2SN')
const formatTime = inject('formatTime')

// ---------------------------------- 暴露方法 ----------------------------------

defineExpose({
  show,
  hide
})
</script>

<style lang="scss" scoped>
.lc-tooltip {
  @apply py-2 px-3 absolute bg-default border rounded z-full;
  box-shadow: rgba(21, 24, 31, 0.16) 0px 12px 24px 0px;
  visibility: visible;
  p {
    @apply text-xs text-default font-semibold;
  }
  .lc-tooltip-item-no-zoom,
  .lc-tooltip-item-zoom {
    @apply flex items-center gap-3;
    &:not(:last-child) {
      @apply mb-1.5;
    }
    .lc-tooltip-color {
      @apply w-5 flex items-center;
    }
    .lc-tooltip-color-rect {
      &::before {
        content: '';
        display: inline-block;
        width: 20px;
        height: 6px;
        border-radius: 2px;
        margin-right: 5px;
        background-color: currentColor;
      }
    }
  }
  .lc-tooltip-tip {
    @apply font-normal text-dimmest text-xs;
  }
}

.lc-tooltip-item-no-zoom {
  .lc-tooltip-step {
    @apply font-semibold;
    &::after {
      content: ':';
      @apply font-semibold;
    }
  }
  .lc-tooltip-value {
    @apply w-10 text-left font-semibold;
  }
  .lc-tooltip-tag {
    @apply truncate;
    max-width: 128px;
  }
}

.lc-tooltip-item-zoom {
  .lc-tooltip-step {
    @apply w-7;
  }
  .lc-tooltip-value {
    @apply col-span-1;
    @apply w-10 text-left;
  }
  .lc-tooltip-time {
    @apply w-28;
  }

  .lc-tooltip-tag {
    @apply truncate;
    max-width: 160px;
  }
}
</style>
