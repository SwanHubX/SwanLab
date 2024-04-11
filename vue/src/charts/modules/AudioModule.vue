<template>
  <!-- 音频容器，完成响应式 -->
  <div class="audios-container" :class="setGrid()" ref="audiosRef">
    <div class="border p-2 rounded" v-for="(audio, index) in audios" :key="audio.title">
      <!-- 多实验时，展示实验标题 -->
      <div class="text-xs flex items-center pb-1" :title="audio.exp" v-if="audio.exp">
        <div class="h-2 w-2 rounded-full shrink-0" :style="{ backgroundColor: audio.color }"></div>
        <p class="pl-1 truncate">{{ audio.exp }}</p>
      </div>
      <div class="audio-container" :ref="(el) => handleAddAudio(index, el)">
        <div class="flex items-center mt-2">
          <PlayButton
            v-model="playingList[index]"
            :color="defaultColor"
            @play="handlePlay(index)"
            @pause="handlePause(index)"
          />
          <!-- 当前时间 -->
          <p class="text-sm ml-1" v-if="audioRef[index]">{{ formatTime(index) }}</p>
          <!-- 文件名称 -->
          <p class="text-sm w-full text-center -ml-10">{{ audios[index].caption }}</p>
          <!-- 下载按钮 -->
          <DownloadButton @click="download(index)" />
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 音频模块，用于播放音频，渲染音频波形图，被AudioChart.vue调用
 * @file: AudioModule.vue
 * @since: 2024-01-30 22:44:56
 **/
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import { ref, inject, reactive, onUnmounted, onMounted } from 'vue'
import PlayButton from '../components/PlayButton.vue'
import { debounce } from '@swanlab-vue/utils/common'
import DownloadButton from '../components/DownloadButton.vue'
const props = defineProps({
  // 接受的音频数据，格式为 [{ audioBuffer: AudioBuffer, title: String, tag: String, caption: String } ]
  audios: {
    type: Object,
    required: true
  }
})
const defaultColor = inject('defaultColor')
// const getColor = inject('getColor')

// 所有音频容器
const audiosRef = ref(null)
// 单个音频容器集合, {index:{audio:channels, dom:HTMLDivElement, ...}}
const audioRef = reactive({})
// 用于控制子播放组件的播放状态
const playingList = ref(props.audios.map(() => false))

// ---------------------------------- 将dom元素和tag对应起来 ----------------------------------
const handleAddAudio = (index, el) => {
  if (audioRef[index]) return // 防止重复添加
  audioRef[index] = {
    channels: sample(props.audios[index].audioBuffer),
    dom: el,
    // 音频源
    ...createAudioContext(props.audios[index].audioBuffer),
    // 当前音频播放的位置,百分比，秒
    offset: 0
  }

  // console.log(audioRef, props.audios)
  // 挂载完成后，开始绘制波形图
  addTaskToBrowserMainThread(() => {
    draw(index)
  })
}

// ---------------------------------- 绘制波形的工具函数 ----------------------------------
/**
 * 抽样
 * @param {AudioBuffer} buffer 音频缓存
 * @returns {array} [positives, negatives] 正数据和负数据
 */
const sample = (buffer) => {
  // 左声道数据
  const leftChannel = buffer.getChannelData(0)
  // 右声道数据
  const rightChannel = buffer.numberOfChannels > 1 ? buffer.getChannelData(1) : null
  // 执行采样，例如每隔10个采样点取一个
  // const sampleRate = buffer.sampleRate
  const step = 10
  const sampledLeftChannel = []
  const sampledRightChannel = []
  for (let i = 0; i < leftChannel.length; i += step) {
    sampledLeftChannel.push(leftChannel[i])
    if (rightChannel) {
      sampledRightChannel.push(rightChannel[i])
    }
  }
  return [sampledLeftChannel, sampledRightChannel]
}

const formatTime = (index) => {
  const duration = audioRef[index].duration
  const offset = audioRef[index].offset
  const now = duration * offset
  const minutes = Math.floor(now / 60)
  const seconds = Math.floor(now % 60)
  // duration格式化，例如 01:30
  const minutes2 = Math.floor(duration / 60)
  const seconds2 = Math.floor(duration % 60)

  return (
    `${minutes}:${seconds < 10 ? '0' + seconds : seconds}` +
    '/' +
    `${minutes2}:${seconds2 < 10 ? '0' + seconds2 : seconds2}`
  )
}

// ---------------------------------- 绘制波形图 ----------------------------------

/**
 * @description: 绘制波形图
 * @param {HTMLDivElement} dom 音频容器,dom 将会往里面添加 canvas
 * @param {AudioBuffer} audioBuffer 音频数据, 基于该数据绘制波形图
 * @param {Number} index 音频数据在 props.audios 中的索引
 * @param {Number} r 比例，用于控制波形图的高度，越大越高
 * @param {Number} a 比例，控制关键帧的高度，越大越高
 */
const draw = (index, r = 2, a = 2 / 3) => {
  const dom = audioRef[index].dom
  // dom.style.width = audiosRef.value.offsetWidth / 2 + 'px'
  const offset = audioRef[index].offset
  const channels = audioRef[index].channels
  const color = defaultColor
  // 将color变淡一些,代表已经播放过的部分
  const color2 = color + '99'
  // console.log('positives', positives, 'negatives', negatives)
  // 如果dom的第一个子元素是canvas，则清空画布
  let canvas, ctx
  const ratio = window.devicePixelRatio || 1
  // 定义画布大小，需要考虑到设备问题
  const width = Math.round(dom.offsetWidth * ratio)
  const height = Math.round(dom.offsetHeight * ratio)
  if (dom.firstChild && dom.firstChild.tagName === 'CANVAS') {
    canvas = dom.firstChild
    // 清空画布
    ctx = canvas.getContext('2d')
    ctx.clearRect(0, 0, canvas.width, canvas.height)
    canvas.width = width
    canvas.height = height
    // console.log('clear, 复用canvas')
  } else {
    // 创建 canvas上下文
    canvas = document.createElement('canvas')
    // 设置画布大小
    canvas.width = width
    canvas.height = height
    ctx = canvas.getContext('2d')
    // 为每个canvas添加一个点击事件，用于控制播放
    canvas.addEventListener('click', (event) => {
      handleClickCanvas(event, index)
    })
    // 挂载 canvas
    canvas.style.width = '100%'
    canvas.style.height = '80%'
    // 添加为第一个子元素
    dom.insertBefore(canvas, dom.firstChild)
  }
  let x = 0
  let y = height / 2
  // let offset = 0
  ctx.lineWidth = 1
  ctx.fillStyle = color2
  // 最大高度为画布高度的0.4倍
  const h = height * 0.4
  //额外延伸5个像素点
  const o = 3
  // 判断长度与画布宽度的关系，一像素对应多少数据
  const step = channels[0].length / width
  const baseHeight = 2
  // 第多少个像素之前是变淡的颜色
  const offsetIndex = Math.min(width - 1, Math.max(0, parseInt(offset * width)))
  // 先从左往右绘制数据
  let num

  for (let k = 0; k < width; k++) {
    // 方便调试
    if (k * step > channels[0].length) throw new Error('k*step > channels[0].length')
    ctx.strokeStyle = k <= offsetIndex ? color2 : color
    // 如果当前处于变淡的区域边界，绘制一条长线段
    if (k === offsetIndex) {
      ctx.beginPath()
      ctx.moveTo(x + k, y - h * a)
      ctx.lineTo(x + k, y + h * a)
      ctx.closePath()
      ctx.stroke()
      ctx.fillStyle = color
    }
    // ctx.beginPath()
    // ctx.moveTo(x + k, y)
    // // 将当前k对应的数据绘制到画布上
    // // 设置高度
    num = channels[0][parseInt(k * step)] * h * r
    num = num > 0 ? num + baseHeight + o : num - baseHeight - o
    const rectHeight = Math.abs(num) + o
    const ry = num > 0 ? y - num : y - o

    // ctx.lineTo(x + k, y - num)
    // ctx.closePath()
    // ctx.stroke()

    // 画矩形
    ctx.fillRect(x + k, ry, 2, rectHeight)
  }

  ctx.fill()
}

// 重新绘制所有波形图，用于前端样式的响应式
const redraw = () => {
  for (let i = 0; i < props.audios.length; i++) {
    draw(i)
  }
}

const debounceRedraw = debounce(redraw, 500)

// ---------------------------------- offset更新 ----------------------------------

function createAudioOffestUpdator(index) {
  const now = Date.now()
  let last = now
  const oldOffset = audioRef[index].offset
  let timer = null
  const update = () => {
    // 如果last和now相差小于15ms
    // console.log('Date.now() - last < 50', Date.now() - last < 50)
    if (Date.now() - last < 50) {
      timer = requestAnimationFrame(update)
      return
    }
    last = Date.now()
    audioRef[index].offset = (Date.now() - now) / audioRef[index].duration / 1000 + oldOffset
    // console.log('audioRef[tag].offset', audioRef[tag].offset, audioRef[tag].duration, oldOffset)
    draw(index)
    // 如果offset大于1，则停止更新
    if (audioRef[index].offset >= 1) {
      // console.log('播放结束')
      destory()
      playingList.value[index] = false
      return
    }
    timer = requestAnimationFrame(update)
  }

  const destory = () => {
    cancelAnimationFrame(timer)
    audioRef[index].updator = null
  }

  return {
    start: () => {
      timer = requestAnimationFrame(update)
    },
    destory
  }
}

// ---------------------------------- 控制某个音频的播放或者暂停 ----------------------------------

/**
 * @description: 创建音频上下文
 */
const createAudioContext = (audioBuffer) => {
  // 假设 audioBuffer 是你的 AudioBuffer 对象
  const audioContext = new (window.AudioContext || window.webkitAudioContext)()
  const source = audioContext.createBufferSource()
  // console.log('audioBuffer', audioBuffer)

  // 将 AudioBuffer 设置为源
  source.buffer = audioBuffer
  const duration = audioBuffer.duration
  // 连接到音频上下文的输出
  source.connect(audioContext.destination)
  return { source, duration }
}

const handlePlay = (index) => {
  // console.log('handlePlay', tag, index, audioRef[tag])
  // 如果已经存在source，且没有暂停，则返回
  if (playingList.value[index]) return console.log('已经在播放了')
  if (audioRef[index].offset >= 1) audioRef[index].offset = 0
  // 如果存在source，且已经暂停，则重新生成一个新的source，并且设置offset为上次暂停的位置
  audioRef[index].source = createAudioContext(props.audios[index].audioBuffer).source
  // 将百分比转换为秒
  const offset = audioRef[index].offset * audioRef[index].duration
  // console.log('offset', offset, audioRef[tag].duration, audioRef[tag].offset)

  audioRef[index].source.start(0, offset)
  playingList.value[index] = true
  // 更新offset
  audioRef[index].updator = createAudioOffestUpdator(index)
  audioRef[index].updator.start()
}

const handlePause = (index) => {
  // console.log('handlePause')
  // console.log(audioRef)
  if (!audioRef[index].source) {
    // console.log('audiosRef.value[tag].source', audiosRef.value[tag].source)
    // 如果不存在source，则返回
    return
  }
  // 暂停音频
  audioRef[index].source.stop()
  playingList.value[index] = false
  audioRef[index].updator?.destory()
}

// ---------------------------------- 点击播放处理 ----------------------------------
const handleClickCanvas = (event, index) => {
  // 更新当前点击水平位置的百分比，重新绘制波形图
  audioRef[index].offset = event.offsetX / event.target.offsetWidth
  draw(index)
  // 如果音频正在播放，删除原来的source，重新生成一个新的source，并且设置offset为点击位置的百分比
  console.log('playingList.value[index]', playingList.value[index])
  if (playingList.value[index]) {
    audioRef[index].source.stop()
    playingList.value[index] = false
    audioRef[index].updator?.destory()
    handlePlay(index, index)
  }
}

// ---------------------------------- 响应式处理 ----------------------------------
// 监听audiosRef的宽度变化，重新绘制波形图
const observer = new ResizeObserver(() => {
  debounceRedraw()
})
onMounted(() => {
  observer.observe(audiosRef.value)
})

// ---------------------------------- 组件卸载，停止播放 ----------------------------------

onUnmounted(() => {
  for (const index in audioRef) {
    try {
      audioRef[index].source?.stop()
      audioRef[index].updator?.destory()
    } catch (_) {
      continue
    }
  }
  observer.disconnect()
})

// ---------------------------------- audioBuffer to blob ----------------------------------
const download = (index) => {
  const blob = props.audios[index].audioBlob
  const filename = props.audios[index].title
  // console.log('blob', blob, filename)
  const a = document.createElement('a')
  a.href = URL.createObjectURL(blob)
  a.download = filename
  a.click()
}

// ---------------------------------- 判断元素个数，如果超过1个，就设置grid-cols-2 ----------------------------------
const setGrid = () => {
  if (props.audios.length > 1) {
    return 'grid-cols-2'
  }
}
</script>

<style lang="scss" scoped>
.audios-container {
  @apply w-full grid gap-4;
  .audio-container {
    @apply h-48 w-full flex flex-col relative text-dimmer;
    /* 如果只有一个元素，沾满两格，否则只占一格 */
  }
}
</style>
