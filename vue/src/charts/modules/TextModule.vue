<template>
  <div class="w-full h-full">
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
      <div class="line" v-for="(text, i) in texts[currentIndex]" :key="text + i" v-show="!skeleton">
        <div class="caption" :title="getCaption(i)">{{ getCaption(i) }}</div>
        <div class="text" :title="text">{{ text }}</div>
        <SLIcon
          icon="zoom"
          class="w-5 h-5 p-1 absolute right-3 cursor-pointer bg-default icon hidden transition-all"
          @click="zoom(text, i)"
        />
      </div>
      <div v-if="skeleton && texts[currentIndex]">
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
      v-model="currentPage"
      :min="pages.minIndex"
      :max="pages.maxIndex"
      :bar-color="color"
      @change="turnPage"
      v-if="texts.length > 1"
    />
    <SLModal max-w="1200" v-model="isZoom">
      <TextDetail :data="current" />
    </SLModal>
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
    type: Array,
    default: () => []
  },
  tag: {
    type: String,
    default: ''
  },
  modal: {
    type: Boolean
  }
})

const emits = defineEmits(['getText'])

const color = inject('colors')[0]
const skeleton = ref(false)

/**
 * 获取 caption
 * @param {*} i
 */
const getCaption = (i) => {
  const line = props.data.list[currentIndex.value]
  // 如果一个 step 只有一个 text
  if (props.texts[currentIndex.value]?.length == 1) {
    return line?.more?.caption || '-'
  }
  // 多个 text，通过行索引得到 caption
  return line?.more[i]?.caption || '-'
}

// ---------------------------------- 分页 ----------------------------------

const pages = computed(() => {
  const minIndex = Math.min(...indexes.value)
  const maxIndex = Math.max(...indexes.value)
  return { minIndex, maxIndex, sum: indexes.value.length }
})

const indexes = computed(() => {
  return props.data.list.map((item) => item.index)
})

// 当前页码，与 step 对应
const currentPage = ref(pages.value.minIndex)
const previousPage = ref(currentPage.value)
// 当前索引，与数据在数组中的索引对应
const currentIndex = ref(0)

const findClosestNumber = (targetNumber, larger) => {
  // 将字符串转换为数字
  const numericArray = indexes.value.map(Number)
  // 计算每个数字与给定数字之间的距离
  const distances = numericArray.map((num) => Math.abs(num - targetNumber))
  // 找到距离最近的数字
  const minDistance = Math.min(...distances)
  let index = distances.indexOf(minDistance)
  let number = numericArray[index]
  if (number === targetNumber) return { index, number }
  // 这个时候已经找到绝对值上最接近的数，但是需要判断是向前还是向后翻页
  // if (
  //   (larger && number <= currentPage.value && index !== indexes.value.length - 1) ||
  //   (!larger && number > currentPage.value && index !== 0)
  // ) {
  //   index += larger ? 1 : -1
  //   return { index, number: numericArray[index] }
  // }

  return { index, number }
}

/**
 * 翻页时展示骨架屏
 * @param {*} tag
 * @param {*} index
 */
const time = ref()
const turnPage = (p) => {
  const { index, number } = findClosestNumber(p, p > previousPage.value)
  currentIndex.value = index
  currentPage.value = Number(number)
  previousPage.value = Number(number)
  skeleton.value = true
  if (time.value) {
    clearTimeout(time.value)
  }
  time.value = setTimeout(() => {
    skeleton.value = false
  }, 400)
  // debounce(() => {
  //   skeleton.value = false
  // }, 300)()
  emits('getText', props.tag, index)
}

// // ---------------------------------- 放大 ----------------------------------

const isZoom = ref(false)

const current = ref({})
/**
 * 放大查看行详情
 * @param {*} text
 * @param {*} i
 */
const zoom = (text, i) => {
  isZoom.value = true
  const line = props.data?.list[currentIndex.value]
  current.value = {
    tag: props.tag,
    line,
    caption: Array.isArray(line?.more) ? line?.more[i]?.caption : line?.more?.caption,
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
