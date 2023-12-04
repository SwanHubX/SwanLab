import { createApp } from 'vue'
import './style.css'
import App from './App.vue'
import './theme/color.min.css'
import { createPinia } from 'pinia'
import { i18n } from './i18n'

/********************************************
 * 颜色主题，从localStorage中读取用户设置的深浅主题偏好，如果没有设置，则默认
 *******************************************/
const theme = localStorage.getItem('theme') || 'dark'
if (theme) {
  // 给body添加data-theme属性
  document.body.setAttribute('data-theme', theme)
} else {
  // 默认白天模式
  document.body.setAttribute('data-theme', 'light')
}

// 設置 pinia
const pinia = createPinia()
const app = createApp(App)

app.use(pinia)
app.use(i18n)
app.mount('#app')
