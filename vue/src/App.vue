<template>
  <EChartsTest :data="data" />
</template>

<script setup>
import EChartsTest from './components/EChartsTest.vue'
import http from './api/http'
import { ref } from 'vue'
// 折线图集合
const data = ref([5, 20, 36])
// 每隔1秒钟，将data中随机的一项改变，增加+5
setInterval(() => {
  http.get('http://localhost:10101/api/test').then((res) => {
    data.value = [...data.value, res.data]
    console.log(data.value)
  })
}, 2000)
</script>

<style scoped>
.logo {
  height: 6em;
  padding: 1.5em;
  will-change: filter;
  transition: filter 300ms;
}
.logo:hover {
  filter: drop-shadow(0 0 2em #646cffaa);
}
.logo.vue:hover {
  filter: drop-shadow(0 0 2em #42b883aa);
}
</style>
