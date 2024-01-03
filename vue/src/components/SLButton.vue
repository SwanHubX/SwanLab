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

// 反转默认样式
.sa-button-reverse {
  background-color: var(--foreground-default);
  border-color: var(--outline-default);
  color: var(--background-default);

  &:hover {
    background-color: var(--foreground-dimmer);
  }

  &:active {
    border-color: var(--primary-default);
  }
}

// 空心默认样式
.sa-button-default-hollow {
  color: var(--foreground-dimmer);
  border-color: var(--outline-default);
  color: var(--foreground-default);
  &:hover {
    background-color: var(--background-hightest);
  }
  &:active {
    border-color: var(--primary-default);
  }
}

// 渐变默认样式
.sa-button-default-gradient {
  color: var(--foreground-default);
  background: linear-gradient(180deg, var(--background-stronger) 0%, var(--background-higher) 100%);
  border-color: transparent;

  &:hover {
    // 取消渐变
    background: var(--background-stronger);
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

// 空心primary样式
.sa-button-primary-hollow {
  color: var(--primary-strongest);
  border-color: var(--primary-dimmer);
  &:hover {
    background-color: var(--primary-dimmest);
  }
  &:active {
    border-color: var(--foreground-default);
  }
}

// 渐变primary样式
.sa-button-primary-gradient {
  background: linear-gradient(180deg, var(--primary-default) 0%, var(--primary-dimmer) 100%);
  border-color: var(--primary-dimmer);
  color: var(--foreground-default);
  &:hover {
    // 取消渐变
    background: var(--primary-default);
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

// 空心positive样式
.sa-button-positive-hollow {
  color: var(--positive-stronger);
  border-color: var(--positive-dimmer);
  &:hover {
    color: var(--positive-strongest);
    background-color: var(--positive-dimmest);
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

// 空心negative样式
.sa-button-negative-hollow {
  color: var(--negative-stronger);
  border-color: var(--negative-dimmer);
  &:hover {
    color: var(--negative-strongest);
    background-color: var(--negative-dimmest);
  }
  &:active {
    border-color: var(--primary-default);
  }
}

// 渐变negative样式
.sa-button-negative-gradient {
  background: linear-gradient(180deg, var(--negative-stronger) 0%, var(--negative-default) 100%);
  color: var(--foreground-default);
  border-color: transparent;
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
