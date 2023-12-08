import { defineConfig, loadEnv } from 'vite'
import { viteMockServe } from 'vite-plugin-mock'
import vue from '@vitejs/plugin-vue'
import path from 'path'
// https://vitejs.dev/config/
export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, path.resolve(process.cwd(), 'vue'))
  // console.log('当前模式：', mode)
  console.log('当前环境：', env)
  const useMock = mode === 'mock'
  return {
    // 服务插件
    plugins: [
      vue(),
      viteMockServe({
        mockPath: 'vue/src/mock',
        localEnabled: useMock
      })
    ],
    root: 'vue',
    // 重定向
    resolve: {
      alias: {
        '@swanlab-vue': path.resolve(__dirname, 'vue/src')
      }
    },
    // 标明编译后存放的位置
    build: {
      outDir: path.resolve(__dirname, 'swanlab/template'),
      emptyOutDir: true
    },
    // 服务配置
    server: {
      proxy: useMock
        ? undefined
        : {
            '/api': {
              target: env.VITE_SERVER_PROXY,
              changeOrigin: true
            }
          },
      host: '0.0.0.0',
      port: 5175,
      open: '.'
    }
  }
})
