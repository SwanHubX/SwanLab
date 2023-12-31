<template>
  <div class="w-full">
    <h1 class="text-xl font-semibold">{{ $t(`common.config-editor.title.${type}`) }}</h1>
    <div class="relative pt-4">
      <h2 class="font-semibold pb-3">{{ $t(`common.config-editor.sub-title.${type}.name`) }}</h2>
      <input type="text" class="input" v-model="info.name" :placeholder="`edit your ${type} name here`" />
      <!-- 提示信息 -->
      <span class="absolute bottom-[-20px] left-0 text-xs text-negative-default">{{ errors.name }}</span>
    </div>
    <div class="relative">
      <h2 class="font-semibold pt-5 pb-3">{{ $t(`common.config-editor.sub-title.${type}.desc`) }}</h2>
      <textarea
        class="input"
        rows="10"
        v-model="info.description"
        :placeholder="`edit your ${type} description here`"
      ></textarea>
      <!-- 提示信息 -->
      <span class="absolute bottom-[-20px] left-0 text-xs text-negative-default">{{ errors.description }}</span>
    </div>
    <div class="flex justify-end pt-10">
      <button
        class="p-2 rounded bg-primary-default text-white transition-all"
        :class="handling ? 'pointer-events-none cursor-not-allowed' : 'hover:rounded-lg active:opacity-70'"
        @click="save"
      >
        {{ $t('common.config-editor.save') }}
      </button>
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 编辑项目/实验信息得弹窗布局
 * @file: EditorWrap.vue
 * @since: 2023-12-31 10:30:01
 **/
import http from '@swanlab-vue/api/http'
import { useProjectStore, useExperimentStroe } from '@swanlab-vue/store'
import { ref } from 'vue'

const projectStore = useProjectStore()
const experimentStore = useExperimentStroe()

const props = defineProps({
  type: String
})

const emits = defineEmits(['hideModal'])

const info = ref({
  name: props.type === 'project' ? projectStore.name : experimentStore.name,
  description: props.type === 'project' ? projectStore.description : experimentStore.experiment.description
})
const errors = ref({
  name: '',
  description: ''
})

// ---------------------------------- 重新设置 ----------------------------------

const handling = ref(false)

const save = async () => {
  if (info.value.name === projectStore.name && info.value.description === projectStore.description) {
    return (errors.value.name = 'nothing changed')
  }

  handling.value = true
  await (props.type === 'project' ? handleProject : handleExperiment)()
  emits('hideModal')
  handling.value = false
}

// 设置项目信息
const handleProject = async () => {
  const { data } = await http.patch('/project/update', info.value)
  projectStore.setProject(data.project)
}

// 设置实验信息
const handleExperiment = async () => {
  const { data } = await http.patch(`/experiment/${experimentStore.id}/update`, info.value)
  experimentStore.setExperiment(data.experiment)
}
</script>

<style lang="scss" scoped>
.input {
  @apply w-full p-2 text-sm outline-none border rounded-lg bg-dimmer;
}
</style>
