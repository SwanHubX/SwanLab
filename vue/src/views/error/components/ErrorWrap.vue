<template>
  <section class="w-full h-full flex flex-col justify-between">
    <h1 class="text-2xl font-semibold pb-6">
      <slot name="title"></slot>
    </h1>
    <div class="pb-6">
      <p>{{ $t('error.wrap.reason-title') }}</p>
      <slot></slot>
    </div>
    <p class="pb-6">
      <span>{{ $t('error.wrap.message[0]') }}</span>
      <a :href="document" target="_blank" class="text-primary-default underline">{{ $t('error.wrap.message[1]') }}</a>
      。
    </p>
    <p class="pb-6">{{ $t('error.wrap.error-code') + error_code }}</p>
    <p>{{ $t('error.wrap.time', { time }) }}</p>
  </section>
</template>

<script setup>
/**
 * @description: 错误信息
 * @file: ErrorWrap.vue
 * @since: 2023-12-15 15:19:55
 **/
import { formatTime } from '@swanlab-vue/utils/time'
import { inject } from 'vue'

const time = formatTime(new Date().getTime() - 8 * 60 * 60 * 1000)
const document = 'https://geektechstudio.feishu.cn/wiki/UInBw9eaziv17IkwfrOcHCZ1nbc'

const error_code = inject('error_code')
</script>

<style lang="scss" scoped></style>
