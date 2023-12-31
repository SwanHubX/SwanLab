<template>
  <div class="w-full">
    <h1 class="text-xl font-semibold">{{ $t(`common.config-editor.title.${type}`) }}</h1>
    <div class="relative pt-4">
      <h2 class="font-semibold pb-3">{{ $t(`common.config-editor.sub-title.${type}.name`) }}</h2>
      <input type="text" class="input" v-model="name" :placeholder="`edit your ${type} name here`" />
      <!-- 提示信息 -->
      <span class="absolute bottom-[-20px] left-0 text-xs text-negative-default">{{ errors.name }}</span>
    </div>
    <div class="relative">
      <h2 class="font-semibold pt-5 pb-3">{{ $t(`common.config-editor.sub-title.${type}.desc`) }}</h2>
      <textarea class="input" rows="10" v-model="desc" :placeholder="`edit your ${type} description here`"></textarea>
      <!-- 提示信息 -->
      <span class="absolute bottom-[-20px] left-0 text-xs text-negative-default">{{ errors.desc }}</span>
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
import { useProjectStore } from '@swanlab-vue/store'
import { ref } from 'vue'

const projectStore = useProjectStore()

const props = defineProps({
  type: String
})

const name = ref(projectStore.name)
const desc = ref(projectStore.description)
const errors = ref({
  name: '',
  desc: ''
})

// ---------------------------------- 重新设置 ----------------------------------

const handling = ref(false)

const save = () => {
  if (name.value === projectStore.name && desc.value === projectStore.description)
    return (errors.value.name = 'nothing changed')

  props.type === 'project' ? handleProject() : handleExperiment()
}

const handleProject = () => {}

const handleExperiment = () => {}
</script>

<style lang="scss" scoped>
.input {
  @apply w-full p-2 text-sm outline-none border rounded-lg bg-dimmer;
}
</style>
