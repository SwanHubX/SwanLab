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
    <div class="w-full h-[310px] overflow-y-auto" :class="{ 'h-[60vh]': modal, 'border-b': data.length > 1 }">
      <!-- line -->
      <div class="line" v-for="(text, i) in tableData[currentIndex]?.data" :key="text + i" v-show="!skeleton">
        <!-- caption -->
        <div class="caption">{{ tableData[currentIndex]?.more[i].caption }}</div>
        <!-- text -->
        <div class="text" :title="text">{{ text }}</div>
        <!-- zoom icon -->
        <SLIcon
          icon="zoom"
          class="w-5 h-5 p-1 absolute right-3 cursor-pointer bg-default icon hidden transition-all"
          @click="zoom(text, i)"
        />
      </div>
      <!-- 骨架屏 -->
      <div v-if="skeleton && data[currentIndex]">
        <div class="flex items-center border-b border-dimmest" v-for="i in [1, 2, 3]" :key="i">
          <div class="md:w-40 w-24 h-10 px-4 shrink-0 flex items-center border-r">
            <span class="skeleton w-full h-1/2"></span>
          </div>
          <div class="flex items-center w-full h-10 ml-4"><span class="skeleton w-1/2 h-1/2"></span></div>
        </div>
      </div>
    </div>
    <!-- 翻页 -->
    <SlideBar
      class="pt-2"
      v-model="currentPage"
      :min="pages.minIndex"
      :max="pages.maxIndex"
      :bar-color="color"
      @change="turnPage"
      @turn="clickToTurn"
      :key="pages.maxIndex"
      v-if="data.length > 1"
      :turn-by-arrow="modal && !isZoom"
    />
    <!-- 数据详情 -->
    <SLModal max-w="1200" v-model="isZoom" esc-exit>
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
import { ref, inject, computed, watch } from 'vue'
import SLModal from '@swanlab-vue/components/SLModal.vue'
import TextDetail from './TextDetail.vue'
import SlideBar from '../components/SlideBar.vue'

const props = defineProps({
  data: {
    type: Object
  },
  tag: {
    type: String,
    default: ''
  },
  modal: {
    type: Boolean
  },
  modelValue: {
    type: Boolean
  }
})

const emits = defineEmits(['getText', 'update:modelValue'])

const color = inject('defaultColor')
const skeleton = ref(false)

const tableData = computed(() => {
  const data = props.data
  return data.map((item) => {
    if (!Array.isArray(item.data)) {
      item.data = [item.data]
      item.more = [item.more]
    }
    return item
  })
})

// ---------------------------------- 分页 ----------------------------------

/**
 * 分页信息，包括 minIndex，maxIndex，sum
 * minIndex: 最小页码
 * maxIndex: 最大页码
 * sum: 总页数
 */
const pages = computed(() => {
  const minIndex = Math.min(...indexes.value)
  const maxIndex = Math.max(...indexes.value)
  return { minIndex, maxIndex, sum: indexes.value.length }
})

/**
 * 对下面三个变量：
 * indexes: 对应 steps
 * currentPage: 当前页码，与 step 对应
 * currentIndex: 当前索引，与数据在数组中的索引对应
 *
 * 这样做的原因是，step 并不一定均匀，而是单增随机的
 * 所以 currentPage 是当前页码
 * 而 currentIndex 是当前页码所对应的数组索引
 * 对应到 props.data 的 list 中，currentIndex 就是数组索引，currentPage 就是索引对应元素中，index 属性值
 *
 * 这二者的关系被抽离，通过 indexes 联系：
 * indexes 中，把 props.data.list 中每一项的 step 提出，按照原本顺序排列
 * 当有 currentIndex 时，可以通过 indexes[currentIndex] 找到数据对应的 step
 * 而有 currentPage 时，也可以通过找到其在 indexes 中的位置而知道 index,从而通过 props.data.list[index] 获取数据
 */
const indexes = computed(() => {
  return props.data.map((item) => item.index)
})
const currentPage = ref(pages.value.minIndex)
const currentIndex = ref(0)

onMounted(() => {
  currentIndex.value = Math.max(indexes.value.length - 1, 0)
  currentPage.value = indexes.value[currentIndex.value]
})

/**
 * 在翻页时找到有效页码和页码索引
 * @param {int} targetNumber 翻页页码
 * @param {boolean} next 当前页码是向前还是向后, true 为向后翻页
 * @param {boolean} isClick 当前翻页操作是否是通过点击上下翻页按钮触发
 */
const findClosestNumber = (targetNumber, next, isClick) => {
  // 将字符串转换为数字
  const numericArray = indexes.value.map(Number)
  // 计算每个数字与给定数字之间的距离
  const distances = numericArray.map((num) => Math.abs(num - targetNumber))
  // 找到距离最近的数字
  const minDistance = Math.min(...distances)
  let index = distances.indexOf(minDistance)
  let number = numericArray[index]
  if (number === targetNumber) return { index, number }
  if (!isClick) return { index, number }
  // 这个时候已经找到绝对值上最接近的数，但是需要判断是向前还是向后翻页
  if ((next && index !== indexes.value.length - 1) || (!next && index !== 0)) {
    index += next ? 1 : -1
    return { index, number: numericArray[index] }
  }

  return { index, number }
}

/**
 * 翻页，同时展示骨架屏
 * @param {number} p 当前页面
 * @param {boolean} isClick 是否是通过点击上下翻页按钮而触发
 */
const time = ref()
const turnPage = (p, isClick) => {
  // 获取有效页码，index 是页码对应数据数组中的索引，number 是页码
  const { index, number } = findClosestNumber(p, isClick, isClick !== undefined)
  currentIndex.value = index
  currentPage.value = Number(number)
  // 骨架屏
  skeleton.value = true
  if (time.value) clearTimeout(time.value)
  time.value = setTimeout(() => {
    skeleton.value = false
  }, 400)
}

/**
 * 通过点击上下按钮翻页
 * @param {string} type 前后，前为 forward，页码增加，后为 backward，页码减小
 * @param {int} p 当前页码
 */
const clickToTurn = (type, p) => {
  if (type === 'forward') return turnPage(p + 1, true)
  return turnPage(p - 1, false)
}

// ---------------------------------- 放大 ----------------------------------

const isZoom = ref(false)

const current = ref({})
/**
 * 放大查看行详情
 * @param {string} text 文本内容
 * @param {int} i v-for 时的行号索引
 *
 * 因为一个 step 下可能有多条文本，而单行和多行时对应数据结构并不一样，所以需要分辨一下
 */
const zoom = (text, i) => {
  // 当前页面所有的信息
  const line = props.data[currentIndex.value]
  current.value = {
    tag: props.tag,
    line,
    // 单行时，more.caption 是 string，多行时为 array
    caption: Array.isArray(line?.more) ? line?.more[i]?.caption : line?.more?.caption,
    text
  }
  isZoom.value = true
}

watch(
  () => isZoom.value,
  (v) => {
    emits('update:modelValue', v)
  }
)
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
