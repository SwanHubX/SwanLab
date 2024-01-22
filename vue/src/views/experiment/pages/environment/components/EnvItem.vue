<template>
  <div class="w-full flex" v-if="envValue">
    <!-- key -->
    <span class="text-dimmer w-44 shrink-0">{{ $t(`experiment.env.keys.${envKey}`) }}</span>
    <span
      class="break-all relative"
      :class="[highLight ? 'high-light' : '', copy ? 'copy' : '']"
      @mouseenter="hover = true"
      @mouseleave="hover = false"
    >
      <span v-if="!link">{{ envValue }}</span>
      <a :href="envValue" target="_blank" class="hover:underline underline-offset-2" v-else>{{ envValue }}</a>
      <SLCopy
        :text="envValue"
        class="absolute top-1 right-1 transition-all opacity-0"
        :class="hover ? 'copy-button' : ''"
        v-show="hover"
      />
    </span>
  </div>
</template>

<script setup>
/**
 * @description: 一行分为两部分：env 和 value，主要是展示作用
 * @file: EnvItem.vue
 * @since: 2024-01-09 19:44:16
 **/

import { ref } from 'vue'

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
    default: true
  }
})

const hover = ref(false)
</script>

<style lang="scss" scoped>
.high-light {
  @apply text-dimmer bg-higher px-2;
}

.copy {
  @apply transition-all overflow-hidden;
  &:hover {
    @apply bg-higher pr-8;
  }
}

.copy-button {
  @apply opacity-100;
}
</style>
