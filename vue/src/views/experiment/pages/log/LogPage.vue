<template>
  <div class="w-full h-full px-7 py-6">
    <section class="log-area">
      <div class="log-line" v-for="i in 1000" :key="i">
        <!-- 行数 -->
        <span>{{ i }}</span>
        <!-- 日志内容 -->
        <span> dsadsadasdsadsa </span>
      </div>
    </section>
  </div>
</template>

<script setup>
/**
 * @description: 实验-日志页
 * @file: LogPage.vue
 * @since: 2023-12-09 20:40:44
 **/
</script>

<style lang="scss" scoped>
.log-area {
  @apply bg-dimmer w-full h-full overflow-auto rounded p-4;
  font-size: 13px;
  line-height: 16px;
  font-family: 'JetBrains Mono', monospace;
  letter-spacing: 0.1px;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  .log-line {
    @apply flex gap-2;
    span {
      @apply block;
    }
    span:first-child {
      @apply w-8 text-right text-gray-400;
    }
  }
}
</style>
