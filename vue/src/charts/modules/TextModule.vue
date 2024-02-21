<template>
  <div class="w-full h-full pb-5">
    <div class="w-full" v-for="(tag, index) in source" :key="tag">
      <!-- header -->
      <div class="w-full flex items-center bg-higher border-y">
        <div class="caption">Caption</div>
        <div class="text">Text</div>
      </div>
      <!-- body -->
      <div class="w-full border-b">
        <!-- line -->
        <div class="line" v-for="line in slice(data[tag].list, index)" :key="line">
          <!-- caption -->
          <div class="caption" :title="line.more?.caption">{{ line.more ? line.more.caption : '-' }}</div>
          <!-- text -->
          <div class="text" :title="line.data">{{ line.data }}</div>
          <!-- zoom -->
          <SLIcon
            icon="zoom"
            class="w-5 h-5 p-1 absolute right-3 cursor-pointer bg-default icon hidden transition-all"
            @click="zoom(tag, line)"
          />
        </div>
      </div>
      <SlideBar
        class="pt-2"
        v-model="currentPage[index]"
        :max="totalPage[index] ? totalPage[index] : 1"
        :min="1"
        :bar-color="color"
        v-if="currentPage != [] && totalPage[index] && totalPage[index] != 1"
      />
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
import { ref, inject, onMounted } from 'vue'
import SLModal from '@swanlab-vue/components/SLModal.vue'
import TextDetail from './TextDetail.vue'
import SlideBar from '../components/SlideBar.vue'

const props = defineProps({
  data: {
    type: Object,
    default: () => {}
  },
  source: {
    type: Array,
    default: () => []
  }
})

const color = inject('colors')[0]

// ---------------------------------- 分页 ----------------------------------

const currentPage = ref([])
const pageSize = ref(10)
const totalPage = ref([])

onMounted(() => {
  props.source.forEach((tag) => {
    currentPage.value.push(1)
    const temp = Math.ceil(props.data[tag].sum / pageSize.value)
    totalPage.value.push(temp)
  })
})

const slice = (data, index) => {
  return data.slice((currentPage.value[index] - 1) * 10, currentPage.value[index] * 10)
}

// ---------------------------------- 放大 ----------------------------------

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
  @apply w-40 py-2 px-4 border-r truncate;
}

.text {
  @apply w-full py-2 px-4 truncate;
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
