import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'
import path from 'path'
// https://vitejs.dev/config/
export default defineConfig({
  plugins: [vue()],
  root: './vue',
  resolve: {
    alias: {
      '@swanlab-vue': path.resolve(__dirname, 'vue/src')
    }
  },
  build: {
    outDir: path.resolve(__dirname, 'swanlab/template'),
    emptyOutDir: true
  },
  server: {
    host: '0.0.0.0',
    port: 5175,
    open: '.'
  }
})
