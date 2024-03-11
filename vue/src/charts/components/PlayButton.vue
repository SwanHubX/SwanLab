<template>
  <SLButton class="p-1 rounded border flex items-center justify-center" :style="{ color }" @click="handlePlay">
    <SLIcon class="h-5 w-5" icon="pause" v-if="modelValue" />
    <SLIcon class="h-5 w-5" icon="play" v-else />
  </SLButton>
</template>

<script setup>
/**
 * @description: 播放按钮
 * @file: PlayButton.vue
 * @since: 2024-01-30 21:20:05
 **/
import { ref } from 'vue'

const props = defineProps({
  color: {
    type: String,
    required: true
  },
  modelValue: {
    type: Boolean,
    required: true
  }
})

const emit = defineEmits(['play', 'pause', 'update:modelValue'])

const handlePlay = () => {
  if (!props.modelValue) {
    emit('play')
    emit('update:modelValue', true)
  } else {
    emit('pause')
    emit('update:modelValue', false)
  }
}
</script>

<style lang="scss" scoped></style>
