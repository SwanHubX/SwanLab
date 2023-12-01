import { defineConfig, loadEnv } from 'vite'
import vue from '@vitejs/plugin-vue'
import path from 'path'
// https://vitejs.dev/config/
export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, process.cwd())
  return {
    plugins: [vue()],
    root: './vue',
    envDir: './vue',
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
      proxy: {
        '/api': {
          target: 'http://127.0.0.1:10101',
          changeOrigin: true,
        },
      },
      host: '0.0.0.0',
      port: 5175,
      open: '.'
    }
  }
})
