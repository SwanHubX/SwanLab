<template>
  <div class="flex-col justify-center items-center pt-10">
    <div class="loader"></div>
    <p class="text-center pt-10 text-lg">{{ $t('home.empty.table') }}</p>
  </div>
</template>

<script setup>
/**
 * @description: 没有实验时的占位符
 * @file: EmptyTable.vue
 * @since: 2024-01-04 20:26:21
 **/
</script>

<style lang="scss" scoped>
.loader {
  width: 106px;
  height: 56px;
  display: block;
  margin: 30px auto;
  background-image: linear-gradient(#333 50px, transparent 0), linear-gradient(#333 50px, transparent 0),
    linear-gradient(#333 50px, transparent 0), linear-gradient(#333 50px, transparent 0),
    radial-gradient(circle 14px, #333 100%, transparent 0);
  background-size: 48px 15px, 15px 35px, 15px 35px, 25px 15px, 28px 28px;
  background-position: 25px 5px, 58px 20px, 25px 17px, 2px 37px, 76px 0px;
  background-repeat: no-repeat;
  position: relative;
  transform: rotate(-45deg);
  box-sizing: border-box;
}

.loader::after,
.loader::before {
  content: '';
  position: absolute;
  width: 56px;
  height: 56px;
  border: 6px solid #333;
  border-radius: 50%;
  left: -45px;
  top: -10px;
  background-repeat: no-repeat;
  background-image: linear-gradient(#333 64px, transparent 0), linear-gradient(#333 66px, transparent 0),
    radial-gradient(circle 4px, #333 100%, transparent 0);
  background-size: 40px 1px, 1px 40px, 8px 8px;
  background-position: center center;
  box-sizing: border-box;
  animation: rotation 0.3s linear infinite;
}

.loader::before {
  left: 25px;
  top: 60px;
}

@keyframes rotation {
  0% {
    transform: rotate(0deg);
  }
  100% {
    transform: rotate(360deg);
  }
}
</style>
