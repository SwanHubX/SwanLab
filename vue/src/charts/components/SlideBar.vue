<template>
  <div class="chart-slide">
    <p>{{ reference }}</p>
    <p>{{ min }}</p>
    <div class="slide">
      <SLSlideBar is-int :max="max" :min="min" v-model="_modelValue" :bar-color="barColor" />
    </div>
    <p>{{ max }}</p>
    <!-- 输入框 -->
  </div>
</template>

<script setup>
/**
 * @description: 封装全局slidebar组件，添加一些其他功能
 * @file: SlideBar.vue
 * @since: 2024-01-30 16:18:31
 **/
import { computed } from 'vue'

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
  }
})

const emits = defineEmits(['update:modelValue', 'change'])

const _modelValue = computed({
  get() {
    return props.modelValue
  },
  set(value) {
    emits('update:modelValue', value)
    emits('change', value)
  }
})
</script>

<style lang="scss" scoped>
.chart-slide {
  @apply flex items-center justify-center flex-wrap w-full gap-2 text-dimmer;
  .slide {
    @apply max-w-[230px] w-full;
  }
}
</style>
