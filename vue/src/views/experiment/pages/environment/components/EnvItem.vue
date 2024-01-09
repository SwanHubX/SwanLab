<template>
  <div class="w-full flex" v-if="envValue">
    <!-- key -->
    <span class="text-dimmer w-44 shrink-0">{{ $t(`experiment.env.keys.${envKey}`) }}</span>
    <span class="break-all" :class="highLight ? 'text-dimmer bg-higher px-2' : ''">{{ envValue }}</span>
  </div>
</template>

<script setup>
/**
 * @description: 一行分为两部分：env 和 value，主要是展示作用
 * @file: EnvItem.vue
 * @since: 2024-01-09 19:44:16
 **/

defineProps({
  envKey: {
    type: String,
    default: ''
  },
  envValue: {
    type: [String, Number],
    default: ''
  },
  highLight: {
    type: Boolean,
    default: false
  }
})
</script>

<style lang="scss" scoped></style>
