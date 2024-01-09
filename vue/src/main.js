import { createApp } from 'vue'
import App from './App.vue'
import { createPinia } from 'pinia'
import { i18n } from './i18n'
import router from './router'
import './style.scss'
// 設置 pinia
const pinia = createPinia()
const app = createApp(App)

app.use(pinia)
app.use(i18n)
app.use(router)

// ---------------------------------- 设置指令 ----------------------------------
// v-tippy
import { directive } from 'vue-tippy'
import 'tippy.js/dist/tippy.css'

app.directive('tippy', directive)

// ---------------------------------- 挂载 ----------------------------------

app.mount('#app')

/********************************************
 * 判断是否是苹果设备，如果不是，默认为windows，设置一些样式
 *******************************************/
import { isApple } from '@swanlab-vue/utils/browser'
if (!isApple) {
  const sidevbarTheme = `
    /* 滚动条整体样式 */
    * ::-webkit-scrollbar {
      width: 8px; /* 设置滚动条宽度 */
    }
    /* 滚动条轨道样式 */
    * ::-webkit-scrollbar-track {
      background-color: var(--background-default);
    }

    /* 滚动条滑块样式 */
    * ::-webkit-scrollbar-thumb {
    background-clip: padding-box;
      background-color: rgba(50, 50, 50, 0.3);
      border-radius: 6px;
      border: 2px solid rgba(0,0,0,0);
    }

    *::-webkit-scrollbar-button {
      display: none;
      height: 0;
      width: 12px
    }

    /*
    * ::-webkit-scrollbar-thumb:hover {
      background-color: var(--outline-stronger);
    }

    * ::-webkit-scrollbar-thumb:active {
      background-color: var(--outline-strongest);
    }
    */
  `
  const style = document.createElement('style')
  style.innerHTML = sidevbarTheme
  document.getElementsByTagName('head')[0].appendChild(style)
}
