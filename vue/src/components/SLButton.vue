<template>
  <button :class="themeClass" :disabled="disabled" :title="disabled ? disabledTip : ''" ref="buttonRef">
    <slot></slot>
    {{ text }}
  </button>
</template>

<script setup>
/**
 * SaButton - 封装一些常用的按钮样式，只对背景、字体颜色和边框颜色做了封装，其他样式需要自己写
 * @prop {string} theme - 按钮主题，可选值: primary、positive、negative、default等，默认值: default
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
  },
  disabled: {
    type: Boolean,
    default: false
  },
  disabledTip: {
    type: String,
    default: 'Not Allowed'
  }
})

const themeClass = computed(() => {
  const theme = `sa-button-${props.theme}${props.hollow ? '-hollow' : ''}`
  return props.disabled ? '' : theme
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
    pointer-events: auto !important;
    cursor: not-allowed !important;
    opacity: 0.5;
  }
}

// 默认样式
.sa-button-default {
  @apply bg-higher border-default text-default;

  &:hover {
    @apply bg-highest;
  }

  &:active {
    @apply border-primary-default;
  }
}

// primary样式
.sa-button-primary {
  @apply bg-primary-higher border-primary-higher text-white-default;

  &:hover {
    @apply bg-primary-default border-primary-default;
  }

  &:active {
    border-color: var(--foreground-default);
  }
}

// positive样式
.sa-button-positive {
  @apply text-white-default border-positive-higher bg-positive-higher;

  &:hover {
    @apply border-positive-default bg-positive-default;
  }

  &:active {
    @apply border-primary-default;
  }
}

// negative样式
.sa-button-negative {
  @apply text-white-default border-negative-higher bg-negative-higher;

  &:hover {
    @apply border-negative-default bg-negative-default;
  }

  &:active {
    @apply border-primary-default;
  }
}

// ----------------------- hollow 样式 -----------------------
.sa-button-negative-hollow {
  @apply text-negative-default bg-default transition-all;

  &:hover {
    @apply text-white-default border-negative-higher bg-negative-higher;
  }

  &:active {
    @apply border-primary-default;
  }
}

.sa-button-primary-hollow {
  @apply text-default bg-default transition-all;

  &:hover {
    @apply text-white-default bg-primary-default border-primary-default;
  }

  &:active {
    @apply opacity-70;
  }
}
</style>
