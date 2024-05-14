<template>
  <SLModal class="pt-5 overflow-hidden" max-w="800" v-model="downloadModal" escExit>
    <p class="text-lg px-5 font-semibold">导出为图片</p>
    <div class="relative p-4">
      <LineChart :index="index" ref="downloadRef" title="123" :chart="chart"></LineChart>
    </div>
    <div class="flex justify-end p-5 pt-0">
      <SLButton theme="primary py-1.5 px-3 rounded">下载图片</SLButton>
    </div>
  </SLModal>
</template>

<script setup>
/**
 * @description:
 * @file: ExportImage.vue
 * @since: 2024-05-14 16:22:09
 **/
import SLModal from '@swanlab-vue/components/SLModal.vue'
import SLButton from '@swanlab-vue/components/SLButton.vue'
import LineChart from '../package/LineChart.vue'

const props = defineProps({
  chart: {
    type: Object,
    required: true
  },
  index: {
    type: Number,
    required: true
  },
  modelValue: {
    type: Boolean,
    default: false
  }
})
const emit = defineEmits(['update:modelValue'])
const data = inject('data')
const downloadRef = ref(null)
const downloadModal = computed({
  get() {
    downloadRef.value?.render(data)
    return props.modelValue
  },
  set(val) {
    emit('update:modelValue', val)
  }
})
</script>

<style lang="scss" scoped></style>
