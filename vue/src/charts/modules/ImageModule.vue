<template>
  <div class="image-detail">
    <div class="text-xs flex items-center justify-center pb-1 w-full" :title="name" v-if="name">
      <div class="h-2 w-2 rounded-full shrink-0" :style="{ backgroundColor: color }"></div>
      <p class="pl-1 truncate max-w-full">{{ name }}</p>
    </div>
    <div class="image-container" @click="$emit('zoom', filename, index)">
      <div class="relative">
        <img :src="imagesData" />
        <DownloadButton class="image-download-button" @click.stop="$emit('download', filename)" />
      </div>
    </div>
    <p class="text-xs w-full text-center truncate" v-if="caption">{{ caption }}</p>
  </div>
</template>

<script setup>
/**
 * @description: 图像模块，负责单个图像显示和暴露下载事件
 * @file: ImageModule.vue
 * @since: 2024-03-17 22:11:10
 **/
import DownloadButton from '../components/DownloadButton.vue'
defineProps({
  filename: String,
  imagesData: String,
  // 图像描述
  caption: String,
  index: [Number, String],
  // 实验名称
  name: String,
  color: String,
  // 是否为多实验
  multi: Boolean
})

defineEmits(['zoom', 'download'])
</script>

<style lang="scss" scoped>
.image-detail {
  @apply relative h-full w-full;
  max-height: 13rem;

  .image-container {
    @apply relative w-full grow inline-flex justify-center cursor-pointer;
    height: calc(100% - 2.5rem);
    img {
      @apply my-auto object-contain;
      height: 100%;
    }
    &:hover .image-download-button {
      @apply block;
    }
  }
}
</style>
