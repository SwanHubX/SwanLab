<template>
  <section class="w-screen h-screen overflow-x-clip">
    <!-- 顶部header -->
    <header class="h-14 w-full">
      <HeaderBar :version="version" />
    </header>
    <!-- 下半部分主要内容，增加一个relative-container的原因是解决按钮在滑动时的隐藏问题 -->
    <div class="relative-container">
      <main class="main-container">
        <!-- 遮罩，淡入淡出 -->
        <transition name="fade" v-if="showSideBar">
          <div class="md:hidden sidebar-overlay" v-if="isSideBarShow" @click="handleClose"></div>
        </transition>
        <!-- 侧边栏 -->
        <div class="sidebar-container bg-default" ref="sidebarRef" v-if="showSideBar">
          <!-- 侧边栏规定宽度 -->
          <div class="w-80 h-full border-r">
            <SideBar />
          </div>
        </div>
        <!-- 右侧主要内容 -->
        <div class="main-content" ref="containerRef">
          <!-- 侧边栏关闭/开启按钮 -->
          <button class="close-button" ref="cbRef" @click="handleClose" v-if="showSideBar">
            <SLIcon icon="sidebar" class="w-full h-full" />
          </button>
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
const props = defineProps({
  version: {
    type: String,
    default: undefined
  },
  showSideBar: {
    type: Boolean,
    default: true
  }
})

// 初始按钮top位置，单位px
const initialTop = 26 + 'px'

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
      if (!props.showSideBar) return
      // 显示
      if (val) {
        sidebarRef.value.style = 'width: 320px;'
        cbRef.value.style.transform = ''
        cbRef.value.style.top = initialTop
      } else {
        sidebarRef.value.style = 'width: 0;border: none;'
        // cb添加transform动画，向左移动233px，旋转180度，向下移动60px
        cbRef.value.style.transform = 'translateX(-260px) rotateY(180deg)'
        handleContainerScroll(false)
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

const handleContainerScroll = (removeAnimation = true) => {
  if (isSideBarShow.value) return
  // 获取滚动距离
  requestAnimationFrame(() => {
    const scrollTop = containerRef.value.scrollTop
    // 按钮top位置为其class中的top值-滚动距离
    cbRef.value.style.top = `calc(26px - ${scrollTop}px)`
    removeAnimation && cbRef.value.classList.remove('close-animation')
  })
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
  @apply absolute z-full outline-none border-none rounded p-1 bg-transparent;
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
