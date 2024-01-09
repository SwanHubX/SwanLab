<template>
  <div class="w-full">
    <div class="px-6 py-4 bg-higher rounded" v-if="dependencies && dependencies.length !== 0">
      <p v-for="line in dependencies" :key="line">{{ line }}</p>
    </div>
    <div class="w-full text-center pt-10" v-else>
      {{ $t('experiment.env.empty.requirements') }}
    </div>
  </div>
</template>

<script setup>
/**
 * @description: 依赖项
 * @file: EnvDependances.vue
 * @since: 2024-01-09 16:01:55
 **/

import { useExperimentStroe } from '@swanlab-vue/store'

const experimentStore = useExperimentStroe()
const dependencies = experimentStore.experiment.system?.dependencies
</script>

<style lang="scss" scoped></style>
