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
        <div class="line" v-for="(text, i) in getTexts(tag, index)" :key="text + i" v-show="!skeleton">
          <!-- caption -->
          <div class="caption" :title="getCaption(tag, index, i)">
            {{ getCaption(tag, index, i) }}
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
        <!-- 骨架 -->
        <div v-if="skeleton || !getTexts(tag, index)">
          <div class="flex items-center border-b border-dimmest" v-for="i in [1, 2, 3]" :key="i">
            <div class="md:w-40 w-24 h-10 px-4 shrink-0 flex items-center border-r">
              <span class="skeleton w-full h-1/2"></span>
            </div>
            <div class="flex items-center w-full h-10 ml-4"><span class="skeleton w-1/2 h-1/2"></span></div>
          </div>
        </div>
      </div>
      <SlideBar
        class="pt-2"
        v-model="currentPage[index]"
        :max="totalPage[index] ? totalPage[index] : 1"
        :min="1"
        :bar-color="color"
        :key="totalPage[index]"
        @change="turnPage(tag, currentPage[index])"
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
import { debounce } from '@swanlab-vue/utils/common'

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

const emits = defineEmits(['getText'])

const color = inject('colors')[0]
const skeleton = ref(false)

const getCaption = (tag, index, i) => {
  const texts = props.texts[tag]
  const line = props.data[tag].list
  const cp = currentPage.value[index] - 1
  if (texts[cp]?.length == 1) {
    return line[cp]?.more?.caption || '-'
  }
  return line[cp]?.more[i]?.caption || '-'
}

const getTexts = (tag, index) => {
  return props.texts[tag][currentPage.value[index] - 1]
}

// ---------------------------------- 分页 ----------------------------------

const currentPage = ref(Array(props.source.length).fill(1))
const pageSize = ref(1)
const totalPage = computed(() => {
  return props.source.map((tag) => {
    return Math.ceil(props.data[tag].sum / pageSize.value)
  })
})

const turnPage = (tag, index) => {
  skeleton.value = true
  debounce(() => {
    skeleton.value = false
  }, 300)()
  emits('getText', tag, index)
}

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
  @apply md:w-40 w-24 h-10 px-4 border-r truncate shrink-0 pt-2;
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

@-webkit-keyframes skeleton-ani {
  0% {
    left: 0;
  }

  to {
    left: 100%;
  }
}

@keyframes skeleton-ani {
  0% {
    left: 0;
  }

  to {
    left: 100%;
  }
}

.skeleton {
  @apply block;
  background-image: linear-gradient(-45deg, #f5f5f5 40%, #fff 55%, #f5f5f5 63%);
  list-style: none;
  background-size: 400% 100%;
  background-position: 100% 50%;
  animation: skeleton-animation 2s ease infinite;
}

@keyframes skeleton-animation {
  0% {
    background-position: 100% 50%;
  }

  100% {
    background-position: 0 50%;
  }
}
</style>
