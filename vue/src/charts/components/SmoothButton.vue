<template>
  <SLMenu class="right-0 w-96" corner-class="left-2.5" not-down>
    <div class="smooth-button">
      <SLIcon class="w-full h-full" icon="smooth" />
    </div>
    <template #pop>
      <div class="flex select-none flex-col w-full px-4 text-sm font-semibold gap-4">
        <div class="flex items-center gap-2">
          <span class="w-20 flex-shrink-0">{{ $t('chart.charts.share.smooth.smoothing') }}</span>
          <SLSlideBar
            show-input
            v-model="_value"
            :max="nowMethod.max"
            :min="nowMethod.min"
            @turn="handleTurn"
            @change="handleChange"
            :disabled="nowMethod.disabled"
            :is-int="nowMethod.isInt"
          />
        </div>
        <div class="flex items-center gap-2">
          <span class="w-20 flex-shrink-0">{{ $t('chart.charts.share.smooth.method') }}</span>
          <SLMenu menu-class="w-full" items-class="overflow-auto bg-default" down down-position="top-2.5 right-2">
            <div class="w-full border p-2 text-left rounded hover:bg-higher">{{ nowMethod.name }}</div>
            <template #pop>
              <SLMenuItems class="method-items">
                <SLMenuItem
                  class="method-item"
                  v-for="method in methods"
                  :key="method.id"
                  @click="handleMethodClick(method)"
                >
                  {{ method.name }}
                </SLMenuItem>
              </SLMenuItems>
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
import SLMenuItems from '@swanlab-vue/components/menu/SLMenuItems.vue'
import SLMenuItem from '@swanlab-vue/components/menu/SLMenuItem.vue'
import { debounce } from '@swanlab-vue/utils/common'

const props = defineProps({
  defaultMethodId: {
    type: Number,
    default: 0
  }
})

// 当method被改变/值被改变时，向外触发事件：smooth
const emit = defineEmits(['smooth'])

const debounceSmooth = debounce(() => {
  emit('smooth', {
    id: nowMethodId.value,
    value: _value.value
  })
}, 300)
// ---------------------------------- 算法配置和切换 ----------------------------------
const methods = [
  {
    id: 0,
    name: 'No Smoothing',
    disabled: true
  },
  {
    id: 1,
    name: 'Time Weighted EMA',
    disabled: false,
    min: 0,
    max: 0.99
  },
  {
    id: 2,
    name: 'Running Average',
    disabled: false,
    min: 1,
    max: 100,
    isInt: true
  },
  {
    id: 3,
    name: 'Gaussian Smoothing',
    disabled: false,
    min: 1,
    max: 100,
    isInt: true
  }
]
// 当前选中的算法
const nowMethodId = ref(props.defaultMethodId)
// 当前选中的算法的配置
const nowMethod = computed({
  get: () => {
    // 找不到就返回第一个
    return methods.find((item) => item.id === nowMethodId.value) || methods[0]
  }
})
// 当前选中算法的值
const _value = ref(0)
// ---------------------------------- 当方法被点击 ----------------------------------
const handleMethodClick = (method) => {
  nowMethodId.value = method.id
  _value.value = method.min || 0
  debounceSmooth()
}
// ---------------------------------- 处理slidebar的change和turn事件 ----------------------------------
// 不需要考虑最大最小值溢出的情况，因为slidebar组件已经处理了
const handleTurn = (direction) => {
  console.log('direction', direction)
  if (direction === 'forward') {
    // 向前+1或+0.01,保留两位小数
    _value.value = nowMethod.value.isInt ? _value.value + 1 : Number((_value.value + 0.01).toFixed(2))
  } else {
    // 向后-1或-0.01
    _value.value = nowMethod.value.isInt ? _value.value - 1 : Number((_value.value - 0.01).toFixed(2))
  }
  debounceSmooth()
}
const handleChange = (value) => {
  _value.value = Number(value.toFixed(2))
  debounceSmooth()
}
</script>

<style lang="scss" scoped>
.smooth-button {
  @apply h-8 w-8 p-1 rounded text-dimmer hover:text-positive-default hover:bg-positive-dimmer;
}

.method-items {
  .method-item {
    &:hover {
      @apply bg-higher;
    }
  }
}
</style>
