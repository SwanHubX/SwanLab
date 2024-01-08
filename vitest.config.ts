import { defineConfig } from 'vitest/config'
import path from 'path'

export default defineConfig({
  test: {
    globals: true,
    alias: {
      '@swanlab-vue': path.resolve(__dirname, 'vue/src')
    }
  }
})
