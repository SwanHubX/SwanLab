<template>
  <Menu as="div" class="relative select-none" :class="props.menuClass" v-slot="{ open, close }">
    <!-- 按钮区域 -->
    <div class="relative w-full items-center flex select-none" ref="button">
      <MenuButton class="w-full" :class="[{ 'sa-popover-button-active': open }]" data-sa-menu>
        <slot :open="open">按钮</slot>
        <SLIcon
          icon="down"
          class="h-4 w-4 absolute transition-transform duration-300"
          :class="[$props.downPosition, { 'rotate-180': open }]"
          v-if="down"
        />
      </MenuButton>
    </div>
    <!-- 弹窗区域 -->
    <div class="absolute z-full" :class="[$props.class]" :style="{ minWidth: menuMinWidth + 'px' }" ref="container">
      <MenuItems
        tabindex=""
        as="div"
        class="sa-menu-pop py-2"
        :class="[$props.itemsClass]"
        :style="{ maxHeight: props.menuMaxHeight + 'px' }"
        data-sa-menu
      >
        <slot name="pop" :open="open" :close="close">弹窗区域</slot>
      </MenuItems>
    </div>
    <!-- 尖尖角 -->
    <div class="sa-menu-corner" :class="[props.cornerClass]" v-if="hasCorner && open" data-sa-menu ref="corner"></div>
  </Menu>
</template>

<script setup>
/**
 * @description: 菜单组件, 用于包裹菜单项
 * @file: SLMenu.vue
 * @since: 2024-03-04 21:30:31
 **/
import { Menu, MenuButton, MenuItems } from '@headlessui/vue'
import SLIcon from '../SLIcon.vue'

const props = defineProps({
  // 是否需要down图标
  down: {
    type: Boolean,
    default: false
  },
  // 右侧down的位置
  downPosition: {
    type: String,
    default: 'right-3 top-3'
  },
  // 菜单样式，影响最顶层的样式
  menuClass: {
    type: String
  },
  // 弹窗的样式，直观感受上，对外层组件使用class的时候应该修改弹出框的样式
  class: {
    type: String,
    default: 'w-full'
  },
  // 弹窗最小宽度
  menuMinWidth: {
    type: String,
    default: '120'
  },
  // 弹窗的最大高度（超出此高度后，需要与itemsClass配合出现滑动条）
  menuMaxHeight: {
    type: String,
    default: '240'
  },
  // 默认颜色，一般用于覆盖
  itemsClass: {
    type: String,
    default: 'bg-default'
  },
  // 尖尖角样式，只要输入值就代表有尖尖角
  cornerClass: {
    type: String,
    default: ''
  },
  // 自动翻转，默认开始
  notFlip: {
    type: Boolean,
    default: false
  }
})

const hasCorner = computed(() => {
  return !!props.cornerClass
})
</script>

<style lang="scss" scoped>
.sa-popover-button-active[data-sa-menu] {
  border-color: var(--primary-default) !important;
}

@mixin corner {
  position: absolute;
  display: block;
  width: 12px;
  height: 12px;
  pointer-events: none;
  content: ' ';
  clip-path: polygon(20% 0, 100% 0, 100% 80%, 80% 100%, 0 20%);
}

.sa-menu-corner[data-sa-menu] {
  @apply z-full border-t border-r -rotate-45;
  @apply bg-default;
  @include corner;
}

.sa-menu-pop[data-sa-menu] {
  @apply border rounded w-full text-sm my-1;
  box-shadow: 0px 4px 10px 0px rgba(0, 0, 0, 0.1);
}
</style>
