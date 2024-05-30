<template>
  <div class="lc-tooltip" ref="toolTipRef" v-show="isShow" :key="key">
    <!-- <div class="lc-tooltip-item-zoom" v-if="detail">
      <p class="lc-tooltip-color font-semibold !text-sm"></p>
      <p class="lc-tooltip-step font-semibold !text-sm">{{ $t('chart.charts.share.step') }}</p>
      <p class="lc-tooltip-value font-semibold !text-sm">{{ $t('chart.charts.share.value') }}</p>
      <p class="lc-tooltip-time font-semibold !text-sm">{{ $t('chart.charts.share.time') }}</p>
      <p class="lc-tooltip-tag font-semibold !text-sm">{{ $t('chart.charts.share.tag') }}</p>
    </div> -->
    <template v-if="detail && items.length">
      <div class="lc-tooltip-item-zoom" v-for="item in items" :key="item.color" :style="{ color: item.color }">
        <!-- 颜色 -->
        <span class="lc-tooltip-color lc-tooltip-color-rect"></span>
        <!-- 步数 -->
        <span class="lc-tooltip-step">{{ item.data.index }}</span>
        <!-- 数据 -->
        <span class="lc-tooltip-value">{{ formatNumber2SN(item.data.data) }}</span>
        <!-- 标签 -->
        <span class="lc-tooltip-tag">{{ item.data.series }}</span>
        <!-- 时间 -->
        <span class="lc-tooltip-time">{{ formatTime(item.data.create_time) }}</span>
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
 * @file: LineChartTooltip.vue
 * @since: 2024-02-24 15:36:59
 **/
import { ref, inject } from 'vue'
import { isApple } from '@swanlab-vue/utils/browser'
import { t } from '@swanlab-vue/i18n'

const props = defineProps({
  // 是否显示详细版
  detail: {
    type: Boolean,
    default: false
  }
})

const tip = isApple ? t('chart.charts.line.copy.apple') : t('chart.charts.line.copy.windows')
// 提示框宽度，用于计算提示框位置，但是不直接使用
const tooltipWidth = props.detail ? 400 : 256
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
  // 去除其中带&smooth后缀的数据，屏蔽展示
  data = data.filter((item) => !item.data.series.includes('&smooth'))
  // console.log('show', data, width, x)
  data.sort((a, b) => b.data.data - a.data.data)
  // key.value = Math.random()
  items.value = data
  isShow.value = true
  const left = parseFloat(x.split('px')[0])
  // console.log(left + tooltipWidth, width)
  if (left + tooltipWidth - 100 > width) {
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
  @apply py-2 px-3 absolute bg-default border rounded z-10;
  min-width: 180px;
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
      @apply w-5 flex items-center flex-shrink-0;
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
    @apply font-normal text-dimmest text-xs flex-shrink-0;
  }
}

.lc-tooltip-item-no-zoom {
  @apply text-xs;
  .lc-tooltip-step {
    @apply font-semibold flex-shrink-0;
    &::after {
      content: ':';
      @apply font-semibold;
    }
  }

  .lc-tooltip-value {
    @apply text-left font-semibold whitespace-nowrap;
  }

  .lc-tooltip-tag {
    @apply truncate;
  }
}

.lc-tooltip-item-zoom {
  // @apply text-xs;
  .lc-tooltip-step {
    @apply w-7;
  }

  .lc-tooltip-value {
    @apply text-left;
  }

  .lc-tooltip-time {
    @apply w-32;
  }

  .lc-tooltip-tag {
    @apply truncate;
  }
}
</style>
