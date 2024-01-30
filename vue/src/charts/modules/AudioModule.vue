<template>
  <!-- 音频容器，完成响应式 -->
  <div class="audios-container" ref="audiosRef">
    <div class="audio-container" :ref="(el) => handleAddAudio(tag, el)" v-for="tag in tags" :key="tag"></div>
  </div>
</template>

<script setup>
/**
 * @description: 音频模块，用于播放音频，渲染音频波形图
 * @file: AudioModule.vue
 * @since: 2024-01-30 22:44:56
 **/
import { ref, computed } from 'vue'

const props = defineProps({
  // 接受的音频数据，格式为 { tag: { data: AudioBuffer, title: String } }
  audios: {
    type: Object,
    required: true
  }
})

const tags = computed(() => Object.keys(props.audios.value))

// 所有音频容器
const audiosRef = ref(null)
// 单个音频容器集合, {tag:{audio:AudioBuffer, dom:HTMLDivElement}}
const audioRef = ref({})

// ---------------------------------- 将dom元素和tag对应起来 ----------------------------------
const handleAddAudio = (tag, el) => {
  if (audioRef.value[tag]) return // 防止重复添加

  audioRef.value[tag] = {
    audio: props.audios.value[tag].data,
    dom: el
  }
  // 挂载完成后，开始绘制波形图
  draw(el, props.audios.value[tag].audio)
}

// ---------------------------------- 绘制波形图 ----------------------------------

/**
 * @description: 绘制波形图
 * @param {HTMLDivElement} dom 音频容器,dom 将会往里面添加 canvas
 * @param {AudioBuffer} audioBuffer 音频数据, 基于该数据绘制波形图
 */
const draw = (dom, audioBuffer) => {
  // 创建 canvas上下文
  const canvas = document.createElement('canvas')
  const canvasCtx = canvas.getContext('2d')

  // 挂载 canvas
  dom.appendChild(canvas)
}

// 重新绘制所有波形图，用于前端样式的响应式
const redraw = () => {
  for (const tag in audioRef.value) {
    draw(audioRef.value[tag].dom, audioRef.value[tag].audio)
  }
}

// ---------------------------------- 响应式 ----------------------------------
</script>

<style lang="scss" scoped></style>
