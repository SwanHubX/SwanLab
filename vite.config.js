import { defineConfig, loadEnv } from 'vite'
import { viteMockServe } from 'vite-plugin-mock'
import vue from '@vitejs/plugin-vue'
import path from 'path'
// https://vitejs.dev/config/
export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, path.resolve(process.cwd(), 'vue'))
  return {
    plugins: [
      vue(),
      viteMockServe({
        mockPath: './vue/src/mock'
      })
    ],
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
          target: env.VITE_SERVER_PROXY,
          changeOrigin: true,
        },
      },
      host: '0.0.0.0',
      port: 5175,
      open: '.'
    }
  }
})
