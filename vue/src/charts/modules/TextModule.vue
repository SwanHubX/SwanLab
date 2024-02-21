<template>
  <div class="w-full h-full">
    <div class="w-full border-b" v-for="tag in source" :key="tag">
      <!-- header -->
      <div class="w-full flex items-center bg-higher border-y">
        <div class="caption">Caption</div>
        <div class="text">Text</div>
      </div>
      <!-- body -->
      <div class="w-full">
        <!-- line -->
        <div class="line" v-for="line in data[tag].list" :key="line">
          <!-- caption -->
          <div class="caption">{{ line.more ? line.more.caption : '-' }}</div>
          <!-- text -->
          <div class="text">{{ line.data }}</div>
          <!-- zoom -->
          <SLIcon
            icon="zoom"
            class="w-5 h-5 p-1 absolute right-3 cursor-pointer icon hidden transition-all"
            @click="zoom(tag, line)"
          />
        </div>
      </div>
      <SLModal class="py-10 overflow-hidden" max-w="1200" v-model="isZoom">
        <TextDetail :data="current" />
      </SLModal>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 文字表格，被TextChart.vue调用
 * @file: TextModule.vue
 * @since: 2024-02-20 20:06:45
 **/
import { ref } from 'vue'
import SLModal from '@swanlab-vue/components/SLModal.vue'
import TextDetail from './TextDetail.vue'

defineProps({
  data: {
    type: Object,
    default: () => {}
  },
  source: {
    type: Array,
    default: () => []
  }
})

const isZoom = ref(false)

const current = ref({})

const zoom = (tag, line) => {
  isZoom.value = true
  current.value = {
    tag,
    line
  }
}
</script>

<style lang="scss" scoped>
.caption {
  @apply w-28 py-2 px-4 border-r;
}

.text {
  @apply w-full py-2 px-4;
}

.line {
  @apply border-b border-dimmest flex items-center relative last:border-none;

  &:hover {
    .icon {
      @apply hover:bg-highest hover:rounded block;
    }
  }
}
</style>
