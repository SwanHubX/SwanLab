<template>
  <section class="w-full h-full flex flex-col justify-between">
    <h1 class="text-2xl font-semibold pb-6">{{ $t('error.init-error.title') }}</h1>
    <ul class="pb-6">
      <li v-for="index in [0, 1, 2]" :key="index">
        {{ $t(`error.init-error.reasons[${index}]`) }}
      </li>
    </ul>
    <p class="pb-6">
      <span>{{ $t('error.init-error.message[0]') }}</span>
      <a :href="document" target="_blank" class="text-primary-default underline">{{
        $t('error.init-error.message[1]')
      }}</a>
      。
    </p>
    <p>{{ $t('error.init-error.time', { time }) }}</p>
  </section>
</template>

<script setup>
/**
 * @description:
 * @file: InitError.vue
 * @since: 2023-12-14 14:49:52
 **/
import { formatTime } from '@swanlab-vue/utils/time'

const time = formatTime(new Date().getTime() - 8 * 60 * 60 * 1000)
const document = 'https://geektechstudio.feishu.cn/wiki/UInBw9eaziv17IkwfrOcHCZ1nbc'
</script>

<style lang="scss" scoped></style>
