<template>
  <button :class="themeClass" ref="buttonRef">
    <slot></slot>
    {{ text }}
  </button>
</template>

<script setup>
/**
 * SaButton - 封装一些常用的按钮样式，只对背景、字体颜色和边框颜色做了封装，其他样式需要自己写
 * @prop {string} theme - 按钮主题，可选值: primary、positive、negative、default、reverse，默认值: default
 * @prop {string} text - 按钮文本，可不填，以插槽形式填入
 *
 *
 * @example
 * <SLButton text="你好" theme="primary" hollow />
 */
import { computed, ref } from 'vue'
const props = defineProps({
  theme: {
    type: String,
    default: 'default'
  },
  hollow: {
    type: Boolean,
    default: false
  },
  gradient: {
    type: Boolean,
    default: false
  },
  text: {
    type: String,
    default: ''
  }
})

const themeClass = computed(() => {
  const theme = `sa-button-${props.theme}`
  return theme
})
const buttonRef = ref(null)
defineExpose({
  button: buttonRef
})
</script>

<style lang="scss" scoped>
button {
  @apply border select-none box-border;
  &:disabled {
    pointer-events: none;
    cursor: not-allowed;
    opacity: 0.5;
  }
}

// 默认样式
.sa-button-default {
  background-color: var(--background-higher);
  border-color: var(--outline-default);
  color: var(--foreground-default);

  &:hover {
    background-color: var(--background-hightest);
  }

  &:active {
    border-color: var(--primary-default);
  }
}

// primary样式
.sa-button-primary {
  background-color: var(--primary-dimmer);
  border-color: var(--primary-dimmer);
  color: var(--accent-white-default);

  &:hover {
    background-color: var(--primary-default);
    border-color: var(--primary-default);
  }

  &:active {
    border-color: var(--foreground-default);
  }
}

// positive样式
.sa-button-positive {
  color: var(--foreground-default);
  border-color: var(--positive-dimmest);
  background-color: var(--positive-dimmest);

  &:hover {
    background-color: var(--positive-dimmer);
  }

  &:active {
    border-color: var(--primary-default);
  }
}

// negative样式
.sa-button-negative {
  color: var(--foreground-default);
  border-color: var(--negative-dimmer);
  background: var(--negative-dimmer);

  &:hover {
    background-color: var(--negative-default);
  }

  &:active {
    border-color: var(--primary-default);
  }
}

// warning样式
.sa-button-warning {
  color: var(--foreground-default);
  border-color: var(--warning-dimmer);
  background: var(--warning-dimmer);

  &:hover {
    background-color: var(--warning-dimmest);
  }

  &:active {
    border-color: var(--primary-default);
  }
}
</style>
