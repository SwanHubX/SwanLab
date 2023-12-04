<template>
  <EChartsTest :data="data" />
  <div
    class="w-full h-10 flex items-center border border-default text-positive-default bg-dimmer border-b-primary-default border-b-4"
  >
    hello world
  </div>
  <div v-for="color in colors" :class="color">{{ color }}</div>
</template>

<script setup>
import EChartsTest from './components/EChartsTest.vue'
import http from './api/http'
import { ref } from 'vue'
// 折线图集合
const data = ref([5, 20, 36])
// 每隔1秒钟，将data中随机的一项改变，增加+5
setInterval(() => {
  http.get('/mock/api/test').then((res) => {
    data.value = [...data.value, res.data]
    console.log(data.value)
  })
}, 2000)

const colors = [
  'text-negative-default',
  'text-negative-dimmest',
  'text-positive-default',
  'text-positive-dimmer',
  'text-positive-higher',
  'text-positive-highest',
  'text-primary-dimmest',
  'text-primary-default',
  'bg-dimmer text-dimmest',
  'border border-default',
  'text-dimmer',
  'text-dimmest'
]
</script>

<style scoped></style>
