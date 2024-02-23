<template>
  <div class="w-full h-full">
    <div class="w-full relative" v-for="(tag, index) in source" :key="tag">
      <!-- title: show in modal -->
      <p class="w-full text-center py-6 text-xl font-semibold" v-if="modal">{{ tag }}</p>
      <!-- header -->
      <div class="w-full flex items-center bg-higher border-y">
        <div class="caption">Caption</div>
        <div class="text">Text</div>
      </div>
      <!-- body -->
      <div class="w-full h-[310px] overflow-y-auto border-b" :class="{ 'h-[60vh]': modal }">
        <!-- line -->
        <div class="line" v-for="(text, i) in texts[tag][currentPage[index] - 1]" :key="text + i">
          <!-- caption -->
          <div
            class="caption"
            :title="
              texts[tag][currentPage[index] - 1]?.length != 1
                ? data[tag].list[currentPage - 1]?.more[i]
                : data[tag].list[currentPage - 1]?.more
            "
          >
            {{
              (texts[tag][currentPage[index] - 1]?.length != 1
                ? data[tag].list[currentPage - 1]?.more[i]
                : data[tag].list[currentPage - 1]?.more) || '-'
            }}
          </div>
          <!-- text -->
          <div class="text" :title="text">{{ text }}</div>
          <!-- zoom -->
          <SLIcon
            icon="zoom"
            class="w-5 h-5 p-1 absolute right-3 cursor-pointer bg-default icon hidden transition-all"
            @click="zoom(tag, i, text)"
          />
        </div>
      </div>
      <SlideBar
        class="pt-2"
        v-model="currentPage[index]"
        :max="totalPage[index] ? totalPage[index] : 1"
        :min="1"
        :bar-color="color"
        :key="totalPage[index]"
        @change="$emit('getText', tag, currentPage[index])"
        v-if="currentPage != [] && totalPage[index] && totalPage[index] != 1"
      />
      <SLModal max-w="1200" v-model="isZoom">
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
import { ref, inject, computed } from 'vue'
import SLModal from '@swanlab-vue/components/SLModal.vue'
import TextDetail from './TextDetail.vue'
import SlideBar from '../components/SlideBar.vue'

const props = defineProps({
  data: {
    type: Object,
    default: () => {}
  },
  texts: {
    type: Object,
    default: () => {}
  },
  source: {
    type: Array,
    default: () => []
  },
  modal: {
    type: Boolean
  }
})

defineEmits(['getText'])

const color = inject('colors')[0]

const getLine = () => {}

// ---------------------------------- 分页 ----------------------------------

const currentPage = ref(Array(props.source.length).fill(1))
const pageSize = ref(1)
const totalPage = computed(() => {
  return props.source.map((tag) => {
    return Math.ceil(props.data[tag].sum / pageSize.value)
  })
})

// ---------------------------------- 放大 ----------------------------------

const isZoom = ref(false)

const current = ref({})

const zoom = (tag, i, text) => {
  isZoom.value = true
  const line = props.data[tag]?.list[currentPage.value - 1]
  current.value = {
    tag,
    line,
    caption: line?.more[i],
    text
  }
}
</script>

<style lang="scss" scoped>
.caption {
  @apply md:w-40 w-24 py-2 px-4 border-r truncate shrink-0;
}

.text {
  @apply w-full py-2 px-4 truncate;
}

.line {
  @apply border-b border-dimmest flex items-center relative;

  &:hover {
    @apply bg-higher;
    .icon {
      @apply border border-dimmest rounded-sm hover:bg-highest hover:rounded block;
    }
  }
}
</style>
