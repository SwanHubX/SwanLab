<template>
  <div class="w-full flex relative" v-if="envValue">
    <!-- key -->
    <span class="text-dimmer w-44 shrink-0">{{ $t(`experiment.env.keys.${envKey}`) }}</span>
    <!-- value -->
    <span class="item" :class="{ 'high-light': highLight, copy: copy }" @click="copy && copyText(envValue)">
      <span v-if="!link">{{ envValue }}</span>
      <a :href="envValue" target="_blank" class="hover:underline underline-offset-2" v-else>{{ envValue }}</a>
      <SLIcon icon="copy" class="copy-button" v-if="copy" />
    </span>
  </div>
</template>

<script setup>
/**
 * @description: 一行分为两部分：env 和 value，主要是展示作用
 * @file: EnvItem.vue
 * @since: 2024-01-09 19:44:16
 **/

import { message } from '@swanlab-vue/components/message'
import { copyTextToClipboard } from '@swanlab-vue/utils/browser'

defineProps({
  envKey: {
    type: [String, Number],
    default: ''
  },
  envValue: {
    type: [String, Number],
    default: ''
  },
  highLight: {
    type: Boolean,
    default: false
  },
  link: {
    type: Boolean,
    default: false
  },
  copy: {
    type: Boolean,
    default: false
  }
})

const copyText = (value) => {
  copyTextToClipboard(value)
  message.success('Copy!')
}
</script>

<style lang="scss" scoped>
.item {
  @apply left-48 break-all;

  &:hover {
    .copy-button {
      @apply opacity-100;
    }
  }

  .copy-button {
    @apply absolute w-4 h-4 top-1.5 right-1 transition-all duration-150 opacity-0;
  }
}

.high-light {
  @apply text-dimmer bg-higher px-2;
}

.copy {
  @apply transition-all;
  &:hover {
    @apply rounded pr-8 pl-2 scale-y-100 bg-highest cursor-pointer;
  }
}
</style>
