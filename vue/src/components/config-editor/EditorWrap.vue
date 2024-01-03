<template>
  <form class="w-full" @submit.prevent="save">
    <h1 class="text-xl font-semibold">{{ $t(`common.config-editor.title.${type}`) }}</h1>
    <!-- 修改名称 -->
    <div class="relative pt-4">
      <h2 class="font-semibold pb-3">{{ $t(`common.config-editor.sub-title.${type}.name`) }}</h2>
      <!-- 名称不可为空，且只能是英文字符、数字、下划线、中划线 -->
      <input
        type="text"
        class="input"
        v-model="info.name"
        :placeholder="`edit your ${type} name here`"
        pattern="[a-zA-Z0-9_\-\u4e00-\u9fa5]*"
        required
      />
      <!-- 提示信息 -->
      <span class="absolute bottom-[-20px] left-0 text-xs text-negative-default">{{ errors.name }}</span>
    </div>
    <!-- 修改描述信息 -->
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
    <!-- 确认按钮 -->
    <div class="flex justify-end pt-10">
      <button
        type="submit"
        class="p-2 rounded bg-primary-default text-white transition-all flex items-center"
        :class="
          handling ? 'pointer-events-none cursor-not-allowed opacity-50 gap-2' : 'hover:rounded-lg active:opacity-70'
        "
      >
        <SLLoading size="4" v-if="handling" />
        <span>{{ $t(`common.config-editor.save.${handling ? 'savimg' : 'default'}`) }}</span>
      </button>
    </div>
  </form>
</template>

<script setup>
/**
 * @description: 编辑项目/实验信息得弹窗布局
 * @file: EditorWrap.vue
 * @since: 2023-12-31 10:30:01
 **/
import { useProjectStore, useExperimentStroe } from '@swanlab-vue/store'
import { ref } from 'vue'
import SLLoading from '../SLLoading.vue'
import { message } from '@swanlab-vue/components/message'
import { t } from '@swanlab-vue/i18n'

const projectStore = useProjectStore()
const experimentStore = useExperimentStroe()

const props = defineProps({
  type: String
})

const emits = defineEmits(['confirm'])

// ---------------------------------- 系统参数 ----------------------------------

// v-model 绑定到两个输入栏
const info = ref({
  name: props.type === 'project' ? projectStore.name : experimentStore.name,
  description: props.type === 'project' ? projectStore.description : experimentStore.description
})

// 错误信息
const errors = ref({
  name: '',
  description: ''
})

// ---------------------------------- 重新设置 ----------------------------------

// 是否在处理中
const handling = ref(false)

/**
 * 保存修改
 *
 * 1. 清除错误信息
 * 2. 判断是否一点东西没改
 * 3. 实验模式时，实验名是不能重复的
 * 4. 触发确认函数并显示处理中
 */
const save = async () => {
  errors.value = {
    name: '',
    description: ''
  }

  if (info.value.name === projectStore.name && info.value.description === projectStore.description) {
    return (errors.value.name = 'nothing changed in project config')
  } else if (info.value.name === experimentStore.name && info.value.description === experimentStore.description) {
    return (errors.value.name = 'nothing changed in experiment config')
  }

  // 实验模式中，校验实验名是否重复
  let duplicated = false
  if (props.type === 'experiment') {
    projectStore.experiments.forEach((expr) => {
      // 如果重复，提示错误信息，但是注意实验名称可以和原来的一样
      if (info.value.name === expr.name && info.value.name !== experimentStore.name) {
        duplicated = true
        errors.value.name = 'The experiment name is duplicated'
      }
    })
  }
  if (duplicated) return

  handling.value = true

  emits('confirm', info.value, () => {
    handling.value = false
    message.success(t('common.config-editor.save.success'))
  })
}
</script>

<style lang="scss" scoped>
.input {
  @apply w-full p-2 text-sm outline-none border rounded-lg bg-higher;
}
</style>
