<template>
  <button class="pannel-button" v-tippy="{ content: tip }" @click="handleIconClick" :disabled="disabled">
    <SLIcon class="h-full aspect-square" :class="pannelIconClass" :icon="pannelIcon" />
  </button>
</template>

<script setup>
import SLIcon from '@swanlab-vue/components/SLIcon.vue'

/**
 * @description: chart的按钮样式
 * @file: PanelButton.vue
 * @since: 2023-12-12 21:43:02
 **/

const props = defineProps({
  icon: {
    type: String,
    required: true
  },
  tip: {
    type: String,
    required: true
  },
  // 点击后的icon名称，可以不传入，这样点击前后的icon不变
  clickIcon: {
    type: String,
    required: false
  },
  // 点击后的icon样式，可以不传入，这样点击前后的icon样式不变
  clickClass: {
    type: String,
    required: false
  }
})
const pannelIcon = ref(props.icon)
const pannelIconClass = ref('')
const disabled = ref(false)

const handleIconClick = () => {
  if (disabled.value) return
  if (props.clickIcon) {
    pannelIcon.value = props.clickIcon
    disabled.value = true
  }
  if (props.clickClass) {
    pannelIconClass.value = props.clickClass
    disabled.value = true
  }
}
</script>

<style lang="scss" scoped>
.pannel-button {
  @apply w-4 h-4;
  @apply text-dimmer h-full p-0.5 rounded hover:bg-positive-dimmest hover:text-positive-highest;
  &:disabled {
    @apply cursor-progress opacity-60;
  }
}
</style>
