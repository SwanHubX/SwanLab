import { defineConfig, loadEnv } from 'vite'
import AutoImport from 'unplugin-auto-import/vite'
import Components from 'unplugin-vue-components/vite'
import json5Plugin from 'vite-plugin-json5'
import vue from '@vitejs/plugin-vue'
import path from 'path'
// https://vitejs.dev/config/
export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, path.resolve(process.cwd(), 'vue'))
  // console.log('当前模式：', mode)
  // console.log('当前环境：', env)
  // 如果使用mock模式，不使用代理
  const proxy = {
    '/api': {
      target: env.VITE_SERVER_PROXY,
      changeOrigin: true
    }
  }
  return {
    // 服务插件
    plugins: [
      json5Plugin(),
      vue(),
      // 自动化导入
      AutoImport({
        imports: ['vue', 'vue-router'],
        dts: 'auto-imports.d.ts',
        eslintrc: {
          enabled: true
        }
      }),
      // 自动导入组件，自定义组件库
      Components({
        // 指定组件所在文件夹的位置
        dirs: ['src/components', 'src/layouts'],
        // 文件扩展名
        extensions: ['vue'],
        // 配置type文件生成位置
        dts: 'components.d.ts'
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
      emptyOutDir: true,
      minify: 'terser',
      // 根据模式应用不同的 terser 配置
      terserOptions: {
        // 生产环境移除console
        compress: {
          drop_console: mode === 'release',
          drop_debugger: mode === 'release'
        }
      }
    },
    // 服务配置
    server: {
      proxy,
      host: '0.0.0.0',
      port: 5175,
      open: '.'
    }
  }
})
