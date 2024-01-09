<template>
  <section class="w-screen h-screen overflow-x-clip">
    <!-- 顶部header -->
    <header class="h-14 w-full">
      <HeaderBar :version="version" />
    </header>
    <!-- 下半部分主要内容，增加一个relative-container的原因是解决按钮在滑动时的隐藏问题 -->
    <div class="relative-container">
      <main class="main-container">
        <!-- 侧边栏关闭/开启按钮 -->
        <button class="close-button" ref="cbRef" @click="handleClose">
          <svg xmlns="http://www.w3.org/2000/svg" class="w-full h-full" viewBox="0 0 16 16" fill="none">
            <path
              d="M6.66661 2.66663C6.48979 2.66663 6.32023 2.73686 6.1952 2.86189C6.07018 2.98691 5.99994 3.15648 5.99994 3.33329C5.99994 3.5101 6.07018 3.67967 6.1952 3.8047C6.32023 3.92972 6.48979 3.99996 6.66661 3.99996H14.6666C14.8434 3.99996 15.013 3.92972 15.138 3.8047C15.263 3.67967 15.3333 3.5101 15.3333 3.33329C15.3333 3.15648 15.263 2.98691 15.138 2.86189C15.013 2.73686 14.8434 2.66663 14.6666 2.66663H6.66661ZM14.6666 7.33329H6.66661C6.48979 7.33329 6.32023 7.40353 6.1952 7.52856C6.07018 7.65358 5.99994 7.82315 5.99994 7.99996C5.99994 8.17677 6.07018 8.34634 6.1952 8.47136C6.32023 8.59639 6.48979 8.66663 6.66661 8.66663H14.6666C14.8434 8.66663 15.013 8.59639 15.138 8.47136C15.263 8.34634 15.3333 8.17677 15.3333 7.99996C15.3333 7.82315 15.263 7.65358 15.138 7.52856C15.013 7.40353 14.8434 7.33329 14.6666 7.33329ZM5.99994 12.6666C5.99994 12.4898 6.07018 12.3202 6.1952 12.1952C6.32023 12.0702 6.48979 12 6.66661 12H14.6666C14.8434 12 15.013 12.0702 15.138 12.1952C15.263 12.3202 15.3333 12.4898 15.3333 12.6666C15.3333 12.8434 15.263 13.013 15.138 13.138C15.013 13.2631 14.8434 13.3333 14.6666 13.3333H6.66661C6.48979 13.3333 6.32023 13.2631 6.1952 13.138C6.07018 13.013 5.99994 12.8434 5.99994 12.6666ZM4.47127 5.80463C4.59271 5.67889 4.65991 5.51049 4.65839 5.33569C4.65687 5.16089 4.58676 4.99369 4.46315 4.87008C4.33955 4.74647 4.17234 4.67636 3.99754 4.67484C3.82274 4.67332 3.65434 4.74052 3.52861 4.86196L0.861939 7.52863C0.736958 7.65364 0.666748 7.82318 0.666748 7.99996C0.666748 8.17674 0.736958 8.34627 0.861939 8.47129L3.52861 11.138C3.65434 11.2594 3.82274 11.3266 3.99754 11.3251C4.17234 11.3236 4.33955 11.2534 4.46315 11.1298C4.58676 11.0062 4.65687 10.839 4.65839 10.6642C4.65991 10.4894 4.59271 10.321 4.47127 10.1953L2.27594 7.99996L4.47127 5.80463Z"
            />
          </svg>
        </button>
        <!-- 遮罩，淡入淡出 -->
        <transition name="fade">
          <div class="md:hidden sidebar-overlay" v-if="isSideBarShow" @click="handleClose"></div>
        </transition>
        <!-- 侧边栏 -->
        <div class="sidebar-container bg-default" ref="sidebarRef">
          <!-- 侧边栏规定宽度 -->
          <div class="w-80 h-full border-r">
            <SideBar />
          </div>
        </div>
        <!-- 右侧主要内容 -->
        <div class="main-content" ref="containerRef">
          <slot></slot>
        </div>
      </main>
    </div>
  </section>
</template>
<script setup>
/**
 * MainContentLayout - 页面布局，包含sidebar和右侧部分，用于展示页面内容
 * 右侧部分可以通过传入的props定义是否显示
 */
import { ref, watch, onMounted, provide, computed } from 'vue'
import HeaderBar from './components/HeaderBar.vue'
import SideBar from './components/SideBar.vue'
import { useRoute } from 'vue-router'
import { onUnmounted } from 'vue'
defineProps({
  version: {
    type: String,
    default: 'unknown'
  }
})

// ---------------------------------- 开启/关闭sidebar ----------------------------------

const cbRef = ref(null)
/**
 * 是否是小屏幕，这意味着isSideBarShow在初始化的时候是否为false，如果为false，那么侧边栏将不会显示
 * 后续考虑改为响应式
 */
const threshold = 768
const isSmallScreen = ref(window.innerWidth < threshold)
window.addEventListener('resize', () => {
  isSmallScreen.value = window.innerWidth < threshold
})

// 是否显示主布局的sideBar，可以监听，主要用于实现MainContent和MainHeader之间的交互动画
const isSideBarShow = ref(!isSmallScreen.value)
// 注入isSideBarShow，内部插槽可以通过inject('isSideBarShow')获取到isSideBarShow的值，但是不可以修改
provide(
  'isSideBarShow',
  computed(() => isSideBarShow.value)
)

// 如果需要修改，通过下面注入的openSideBar和closeSideBar方法完成，这两个方法如果被内部组件调用，自动检测是否是小屏幕，如果不是小屏幕，不执行任何操作
provide('openSideBar', () => {
  if (isSmallScreen.value) {
    isSideBarShow.value = true
  }
})
provide('closeSideBar', () => {
  if (isSmallScreen.value) {
    isSideBarShow.value = false
  }
})

// 关闭按钮点击事件
const handleClose = () => {
  cbRef.value.classList.add('close-animation')
  isSideBarShow.value = !isSideBarShow.value
}

// ---------------------------------- 其他对象 ----------------------------------

const sidebarRef = ref(null)
const containerRef = ref(null)

// ---------------------------------- 注册挂载钩子 ----------------------------------

onMounted(() => {
  watch(
    isSideBarShow,
    (val) => {
      // 显示
      if (val) {
        sidebarRef.value.style = 'width: 320px;'
        cbRef.value.style.transform = ''
        cbRef.value.style.top = '26px'
        // 如果不是小屏幕，为container设置属性
        // if (!isSmallScreen.value) {
        //   containerRef.value.style = ''
        // }
      } else {
        sidebarRef.value.style = 'width: 0;border: none;'
        // cb添加transform动画，向左移动233px，旋转180度，向下移动60px
        cbRef.value.style.transform = 'translateX(-260px) rotateY(180deg)'
        // 获取滚动距离
        const scrollTop = Math.max(containerRef.value.scrollTop, 0)
        // 按钮top位置为其class中的top值-滚动距离
        cbRef.value.style.top = `calc(26px - ${scrollTop}px)`

        // 如果不是小屏幕，为container设置属性
        // if (!isSmallScreen.value) {
        //   containerRef.value.style = 'width: 100vw;'
        // }
      }
    },
    {
      immediate: true
    }
  )
})

// ---------------------------------- 路由更改时，判断是否是小屏幕，如果是，隐藏sidebar ----------------------------------

const route = useRoute()
watch(
  computed(() => route.fullPath),
  () => {
    if (isSmallScreen.value) {
      isSideBarShow.value = false
    }
  }
)

// ---------------------------------- 按钮位置修改，监听下滑距离 ----------------------------------

onMounted(() => {
  containerRef.value.addEventListener('scroll', handleContainerScroll)
  onUnmounted(() => {
    containerRef.value.removeEventListener('scroll', handleContainerScroll)
  })
})

const handleContainerScroll = () => {
  if (isSideBarShow.value) return
  // 获取滚动距离
  const scrollTop = containerRef.value.scrollTop
  // 按钮top位置为其class中的top值-滚动距离
  cbRef.value.style.top = `calc(26px - ${scrollTop}px)`
  cbRef.value.classList.remove('close-animation')
}

// ---------------------------------- 暴露对象 ----------------------------------

defineExpose({
  containerRef
})
</script>

<style lang="scss" scoped>
// 动画持续时间
$duration: 400ms;

// 关闭按钮
$close-button-size: 24px;

// 侧边栏宽度，需要与上面js中设置的宽度一致
$sidebar-width: 320px;

// 除去header的高度
$main-content-height: calc(100vh - 56px);

.relative-container {
  @apply relative overflow-y-hidden;
  height: $main-content-height;
}
.main-container {
  @apply w-screen flex  overflow-auto;
  height: $main-content-height;

  .main-content {
    @apply h-full w-full overflow-y-auto overflow-x-hidden;
    // 宽度变化
    transition: width $duration ease-in-out;
  }
}

.close-button {
  @apply absolute z-full outline-none border-none rounded p-1;
  width: $close-button-size;
  height: $close-button-size;
  left: calc($sidebar-width - $close-button-size - 16px);
  &:hover {
    @apply bg-highest;
  }
}
.close-animation {
  transition-property: translateX translateY;
  transition-duration: $duration;
  transition-timing-function: ease-in-out;
}

// 定义一个动画，$duration秒之内不透明度从1到0再到1
@keyframes closeOpen {
  0% {
    opacity: 1;
  }
  50% {
    opacity: 0;
  }
  100% {
    opacity: 1;
  }
}

.close-button-enter-active,
.close-button-leave-active {
  animation: closeOpen $duration ease-in-out;
}

// 遮罩
.sidebar-overlay {
  @apply absolute w-full  overflow-x-hidden z-50;
  height: $main-content-height;
  background-color: var(--background-overlay);
  animation: fadeIn $duration ease;
}
// sidebar容器
.sidebar-container {
  @apply z-50 box-content  md:h-full overflow-y-auto overflow-x-clip flex-shrink-0 md:static absolute;
  height: $main-content-height;
  // 宽度变化
  transition: width $duration ease-in-out;
  // 将滚动条隐藏
  &::-webkit-scrollbar {
    display: none; /* 隐藏WebKit浏览器的滚动条 */
  }
}

// 淡出
.fade-leave-active {
  transition: opacity 0.5s;
}
.fade-leave-to {
  opacity: 0;
}

// 淡入
@keyframes fadeIn {
  from {
    opacity: 0; /* 从完全透明开始 */
  }
  to {
    opacity: 1; /* 渐变到完全不透明 */
  }
}
</style>
