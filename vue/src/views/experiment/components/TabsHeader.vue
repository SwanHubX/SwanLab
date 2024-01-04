<template>
  <div class="w-full pt-2 pl-6 bg-higher flex border-b h-full">
    <router-link :to="tab" v-for="tab in tabs" :key="tab" class="link" active-class="link-active">
      {{ $t(`experiment.tabs.${tab}`) }}
    </router-link>
  </div>
</template>

<script setup>
/**
 * @description: 实验页导航栏
 * @file: TabsHeader.vue
 * @since: 2023-12-09 21:44:40
 **/
const tabs = ['index', 'chart', 'log']
</script>

<style lang="scss" scoped>
.link {
  @apply block py-4 px-8 mx-1 rounded-t-lg border-t border-x border-transparent truncate;
  margin-bottom: -1px;
  z-index: 1;
  position: relative;
  overflow: hidden;
  transition: color 0.2s ease-in-out;

  &::before {
    content: '';
    z-index: -1;
    position: absolute;
    top: 50%;
    left: 50%;
    width: 1em;
    height: 1em;
    border-radius: 50%;
    transform-origin: center;
    transform: translate3d(-50%, -50%, 0) scale3d(0, 0, 0);
    transition: transform 0.25s ease-in-out;
    background-color: #fdfeff;
  }

  &:hover::before {
    transform: translate3d(-50%, -50%, 0) scale3d(15, 15, 15);
  }
}

.link-active {
  transition: none;
  @apply bg-default pointer-events-none text-positive-higher;
  border-color: inherit !important;
}
</style>
