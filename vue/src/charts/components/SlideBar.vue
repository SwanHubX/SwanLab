<template>
  <div class="chart-slide">
    <p class="w-9">{{ reference }}</p>
    <!-- <p>{{ min }}</p> -->
    <div class="slide">
      <SLSlideBar is-int :max="max" :min="min" v-model="_modelValue" :bar-color="barColor" />
    </div>
    <!-- <p>{{ max }}</p> -->
    <!-- 输入框 -->
    <div class="flex items-center relative">
      <input type="number" :value="_modelValue" ref="input" @keyup.enter="handleChange" @blur="handleChange" />
      <div class="w-3 flex-shrink-0 flex-col flex absolute right-1">
        <SLIcon icon="down" class="w-full h-3 -rotate-180 -mb-1" @click="handleClickUp" />
        <SLIcon icon="down" class="w-full aspect-square" @click="handleClickDown" />
      </div>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 封装全局slidebar组件，添加一些其他功能
 * @file: SlideBar.vue
 * @since: 2024-01-30 16:18:31
 **/
import { computed, ref, onUnmounted } from 'vue'
import SLSlideBar from '@swanlab-vue/components/SLSlideBar.vue'

const props = defineProps({
  max: {
    type: Number,
    required: true
  },
  min: {
    type: Number,
    required: true
  },
  modelValue: {
    type: Number,
    required: true
  },
  barColor: {
    type: String,
    required: true
  },
  reference: {
    type: String,
    default: 'Step'
  },
  // 通过方向键切换
  turnByArrow: {
    type: Boolean,
    default: false
  },
  // 方向键的类型，false：左右键，true：上下键
  verticalArrow: {
    type: Boolean,
    default: false
  }
})

const emits = defineEmits(['update:modelValue', 'change', 'turn'])

const input = ref(null)

const _modelValue = computed({
  get() {
    return props.modelValue
  },
  set(value) {
    if (value < props.min) {
      return (_modelValue.value = props.min)
    } else if (value > props.max) {
      return (_modelValue.value = props.max)
    }
    emits('update:modelValue', value)
    emits('change', value)
    input.value.value = _modelValue.value
  }
})

// ---------------------------------- 上下键增加/减少数字 ----------------------------------

/**
 * 向后渐少
 */
const handleClickDown = () => {
  if (_modelValue.value > props.min) {
    emits('turn', 'backward', _modelValue.value)
  }
}

/**
 * 向前增加
 */
const handleClickUp = () => {
  if (_modelValue.value < props.max) {
    emits('turn', 'forward', _modelValue.value)
  }
}

// ---------------------------------- 当input输入结束时，再次赋值 ----------------------------------

const handleChange = (e) => {
  if (e.target.value == _modelValue.value) return
  _modelValue.value = e.target.value
}

// ---------------------------------- 键盘左右切换 ----------------------------------

const handleKeydown = (e) => {
  if (!props.turnByArrow) return
  if (e.key == (props.verticalArrow ? 'ArrowDown' : 'ArrowLeft')) {
    handleClickDown()
  } else if (e.key == (props.verticalArrow ? 'ArrowUp' : 'ArrowRight')) {
    handleClickUp()
  }
}

if (props.turnByArrow) {
  document.addEventListener('keydown', handleKeydown)
}

onUnmounted(() => {
  if (props.turnByArrow) {
    document.removeEventListener('keydown', handleKeydown)
  }
})
</script>

<style lang="scss" scoped>
.chart-slide {
  @apply flex items-center justify-center w-full gap-2 text-dimmer select-none flex-nowrap;
  .slide {
    @apply max-w-[230px] w-full;
  }
  input[type='number']::-webkit-inner-spin-button,
  input[type='number']::-webkit-outer-spin-button {
    -webkit-appearance: none;
    margin: 0;
  }
  input[type='number'] {
    -moz-appearance: textfield;
    appearance: textfield;
  }
}

input {
  @apply w-16 h-6 pl-1 pr-5 rounded border outline-none bg-transparent text-xs focus:border-primary-default;
}
</style>
