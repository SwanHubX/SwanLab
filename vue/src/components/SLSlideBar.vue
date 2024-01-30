<template>
  <div class="sliding-range" ref="slidingRange">
    <!-- 滑动条容器 -->
    <div class="w-full relative h-1 flex">
      <!-- 滑动条进度 -->
      <div class="sliding-bar" :style="{ backgroundColor: barColor }" ref="slidingBar"></div>
      <!-- 滑动条背景容器 -->
      <div class="sliding-container"></div>
      <!-- 滑动按钮 -->
      <div class="sliding-button" ref="slidingButton"></div>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 滑块组件
 * @file: SLSlideBar.vue
 * @since: 2024-01-30 16:06:28
 **/
import { computed, onMounted, ref, watch } from 'vue'
/**
 * 滑块封装组件
 */
const props = defineProps({
  // 滑块的最小值
  min: {
    type: Number,
    default: 0
  },
  // 滑块的最大值
  max: {
    type: Number,
    default: 100
  },
  // 滑块的当前值，通过v-model双向绑定
  modelValue: {
    default: 0
  },
  // 是否为整数
  isInt: {
    type: Boolean,
    default: false
  },
  // 已经滑动的进度条颜色
  barColor: {
    type: String,
    default: 'var(--primary-default)'
  }
})

watch(
  () => props.modelValue,
  () => {
    // 访问一次slidingValue，触发get
    slidingValue.value
  }
)

const emits = defineEmits(['on-change', 'update:modelValue'])
const slidingBar = ref(null)
const slidingButton = ref(null)
const slidingRange = ref(null)
const slidingValue = computed({
  get() {
    const v = formatValue(props.modelValue)
    setButtonAndBarPosition(v)
    return v
  },
  set(val) {
    if (val === slidingValue.value) return
    if (!val) val = 0
    // console.log(val)
    // 设置滑动按钮的位置
    setButtonAndBarPosition(val)
    // 触发事件
    emits('update:modelValue', val)
  }
})

onMounted(() => {
  setButtonAndBarPosition(slidingValue.value)
  // 监听滑动条的点击、拖动事件
  slidingRange.value.addEventListener('mousedown', (e) => {
    const startX = e.clientX
    const rect = slidingRange.value.getBoundingClientRect()
    const startLeft = rect.left
    // 相减得到偏移量
    const offsetLeft = startX - startLeft
    // 总长度
    const slidingRangeWidth = slidingRange.value.offsetWidth
    // console.log(startX, startLeft)
    // 设置滑动按钮的位置
    // console.log(offsetLeft / slidingRangeWidth)
    // 双向绑定
    const newValue = formatValue((offsetLeft / slidingRangeWidth) * props.max)
    slidingValue.value = newValue
    // 触发事件
    emits('on-change', newValue)

    // 滑动事件
    const mousemove = (e) => {
      const moveX = e.clientX
      let left = moveX - startX + offsetLeft
      const process = left / slidingRangeWidth
      const newValue = formatValue(process * props.max)
      slidingValue.value = newValue
      // 触发事件
      emits('on-change', newValue)
    }
    const mouseup = () => {
      // console.log('mouseup')
      document.removeEventListener('mousemove', mousemove)
      document.removeEventListener('mouseup', mouseup)
    }
    document.addEventListener('mousemove', mousemove)
    document.addEventListener('mouseup', mouseup)
  })
})

// 格式化值
function formatValue(value) {
  if (props.isInt) value = Math.floor(value)
  // 值在最小值和最大值之间
  value = Math.min(Math.max(value, props.min), props.max)
  // 保留两位小数，且如果是整数，不显示小数点‘
  value = Number(value.toFixed(2).replace('.00', ''))
  return value
}

// 设置按钮位置和进度条位置，传入的是value
function setButtonAndBarPosition(value) {
  const process = value / props.max
  // 按钮位置
  const percent = Math.min(process * 100, 100) + '%'
  slidingButton.value.style.left = percent
  // 进度条长度
  slidingBar.value.style.width = percent
}
</script>

<style lang="scss" scoped>
.sliding-range {
  @apply hover:cursor-grab active:cursor-grabbing mx-0.5 py-3;
  // 滑动条进度
  .sliding-bar {
    @apply h-full rounded-l-full;
  }
  // 滑动条背景容器
  .sliding-container {
    @apply h-full rounded-r-full grow;
    background-color: var(--outline-default);
  }
  // 滑动按钮
  .sliding-button {
    @apply h-3 w-3 border-2 border-dimmer rounded-full bg-default top-0 -translate-y-[30%] -translate-x-[50%] absolute;
  }

  // 当滑动条被点击时，滑动按钮的样式
  &:active {
    .sliding-button {
      @apply border-dimmer;
    }
  }
}
</style>
