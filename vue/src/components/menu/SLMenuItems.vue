<template>
  <div class="sa-menu-items">
    <slot></slot>
  </div>
</template>

<script setup></script>

<style lang="scss" scoped>
.sa-menu-items {
  @apply w-full;
  &:not(:last-child) {
    @apply border-b pb-2;
  }
}
</style>
