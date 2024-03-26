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
        :maxlength="max_name_len"
        pattern="^(?!\s+$).*$"
        required
      />
      <!-- 数量限制 -->
      <p
        class="absolute bottom-[-20px] text-xs right-0"
        :class="info.name.length === 0 || info.name.length === max_name_len ? 'text-negative-default' : ''"
      >
        {{ `${info.name.length} / ${max_name_len}` }}
      </p>
      <!-- 提示信息 -->
      <span class="tip text-negative-default">{{ errors.name }}</span>
    </div>
    <!-- 修改描述信息 -->
    <div class="relative">
      <h2 class="font-semibold pt-5 pb-3">{{ $t(`common.config-editor.sub-title.${type}.desc`) }}</h2>
      <textarea
        class="input"
        rows="10"
        v-model="info.description"
        :placeholder="`edit your ${type} description here`"
        :maxlength="max_description_len"
      ></textarea>
      <!-- 数量限制 -->
      <p
        class="absolute bottom-[-20px] text-xs right-0"
        :class="info.description?.length === max_description_len ? 'text-negative-default' : ''"
      >
        {{ `${info.description ? info.description.length : 0} / ${max_description_len}` }}
      </p>
      <!-- 提示信息 -->
      <span class="tip text-negative-default">{{ errors.description }}</span>
    </div>
    <!-- 确认按钮 -->
    <div class="flex justify-end pt-10">
      <SLButton type="submit" theme="primary" class="p-2 rounded flex items-center" :disabled="handling">
        <SLLoading size="4" v-if="handling" />
        <span>{{ $t(`common.config-editor.save.${handling ? 'savimg' : 'default'}`) }}</span>
      </SLButton>
    </div>
  </form>
</template>

<script setup>
/**
 * @description: 编辑项目/实验信息得弹窗布局
 * @file: EditorWrap.vue
 * @since: 2023-12-31 10:30:01
 **/
import { useProjectStore, useExperimentStore } from '@swanlab-vue/store'
import { ref } from 'vue'
import SLLoading from '@swanlab-vue/components/SLLoading.vue'
import { message } from '@swanlab-vue/components/message'
import { t } from '@swanlab-vue/i18n'
import { computed } from 'vue'

const projectStore = useProjectStore()
const experimentStore = useExperimentStore()

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

// ---------------------------------- 参数限制 ----------------------------------

const max_name_len = 100
const max_description_len = 255

const checkName = computed(() => {
  const name = info.value.name
  // 判断字符串长度是否超过100个字符
  if (name.length > max_name_len) return false
  return true
})

/**
 * 描述校验
 * 描述不能超过255个字符，可以包含任何字符，前后空格自动去除
 */
const checkDescription = computed(() => {
  if (!info.value.description) return true
  if (info.value.description?.length > max_description_len) {
    return false
  }
  return true
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
  // 清空错误信息
  errors.value = {
    name: '',
    description: ''
  }

  // 检查格式要求
  if (props.type === 'project' && !checkName.value) {
    message.warning('invalid project name')
    return (errors.value.name = 'too long or invalid characters')
  } else if (!checkName.value) {
    message.warning('invalid experiment name')
    return (errors.value.name = 'too long or invalid characters')
  }
  if (!checkDescription.value) {
    message.warning('invalid description')
    return (errors.value.description = 'invalid description')
  }

  // 判断是否一点没变
  if (info.value.name.trim() === projectStore.name && info.value.description === projectStore.description) {
    return (errors.value.name = 'nothing changed in project config')
  } else if (
    info.value.name.trim() === experimentStore.name &&
    info.value.description === experimentStore.description
  ) {
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

.tip {
  @apply absolute bottom-[-20px] left-0 text-xs;
}

button[type='submit'] {
  &:disabled {
    @apply pointer-events-none cursor-not-allowed opacity-50 gap-2;
  }
}
</style>
