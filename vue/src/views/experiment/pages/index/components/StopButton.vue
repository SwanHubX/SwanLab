<template>
  <button
    class="p-2 rounded-full border-2 transition-all duration-200 ml-2"
    :class="showColor ? 'border-red-500' : ' cursor-not-allowed'"
    @mouseover="() => (hover = true)"
    @mouseout="() => (hover = false)"
    @click="() => (showWarning = true)"
  >
    <div class="w-4 h-4 bg-gray-200 transition-all duration-200" :class="showColor ? 'bg-red-500' : ''"></div>
  </button>
  <SLModal class="px-10 pt-10 pb-5 overflow-hidden" maxW="500" v-model="showWarning">
    <div class="w-full">
      <h1 class="text-lg font-semibold text-negative-default">注意!</h1>
      <p class="py-4">停止后不可恢复，且可能导致不可预知的后果，是否停止该实验？</p>
      <div class="flex justify-end">
        <button class="px-2 py-1 border transition-all rounded-lg hover:bg-red-500 hover:text-white" @click="confirm">
          Confirm
        </button>
      </div>
    </div>
  </SLModal>
</template>

<script setup>
/**
 * @description: 停止实验按钮
 * @file: StopButton.vue
 * @since: 2023-12-30 20:30:30
 **/
import { useExperimentStroe } from '@swanlab-vue/store'
import { ref } from 'vue'
import { computed } from 'vue'
import SLModal from '@swanlab-vue/components/SLModal.vue'

// ---------------------------------- 弹窗相关 ----------------------------------

const experiment = useExperimentStroe().experiment
const status = useExperimentStroe().status

const hover = ref(false) // 是否hover
// 展示hover时可点击样式
const showColor = computed(() => {
  return status === 0 && hover.value
})

const showWarning = ref(false) // 弹窗状态

// ---------------------------------- 确认删除 ----------------------------------

const confirm = () => {
  console.log(experiment)
  showWarning.value = false
}
</script>

<style lang="scss" scoped></style>
