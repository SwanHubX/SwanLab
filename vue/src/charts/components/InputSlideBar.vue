<template>
  <SLMenu corner-class="left-3" menu-class="ml-2" menu-min-width="200">
    <template #default="{ open }">
      <input
        type="number"
        class="w-16 text-sm"
        :value="value"
        @click="(e) => handleClickInput(e, open)"
        @change="handleInputChange"
        @keydown.enter="handleInputChange"
      />
    </template>
    <template #pop>
      <div class="px-4">
        <SLSlideBar is-int v-model="value" :max="Number(max)" :min="Number(min)" />
      </div>
    </template>
  </SLMenu>
</template>

<script setup>
import SLSlideBar from '@swanlab-vue/components/SLSlideBar.vue'
import SLMenu from '@swanlab-vue/components/menu/SLMenu.vue'

/**
 * @description: 数字输入滑动条
 * @file: InputSlideBar.vue
 * @since: 2024-05-17 18:06:08
 **/

const props = defineProps({
  modelValue: {
    type: Number,
    required: true
  },
  min: {
    type: String,
    required: true
  },
  max: {
    type: String,
    required: true
  }
})

const emit = defineEmits(['update:modelValue'])

const value = computed({
  get: () => props.modelValue,
  set: (val) => {
    if (val < Number(props.min)) {
      emit('update:modelValue', Number(props.min))
    } else if (val > Number(props.max)) {
      emit('update:modelValue', Number(props.max))
    } else {
      emit('update:modelValue', val)
    }
  }
})

// ---------------------------------- handleInputChange ----------------------------------
const handleInputChange = (e) => {
  const val = Number(e.target.value)
  if (val < Number(props.min)) {
    emit('update:modelValue', Number(props.min))
  } else if (val > Number(props.max)) {
    emit('update:modelValue', Number(props.max))
  } else {
    emit('update:modelValue', val)
  }
}

// ---------------------------------- handleClickInput ----------------------------------
const handleClickInput = (e, open) => {
  // console.log('handleClickInput', e, open)
  if (open) {
    // console.log('阻止关闭')
    e.stopPropagation()
  }
}
</script>

<style lang="scss" scoped>
input {
  @apply border px-2 py-1 rounded;
}
</style>
