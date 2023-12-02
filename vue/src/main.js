import { createApp } from 'vue'
import './style.css'
import App from './App.vue'
import './theme/color.min.css'

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

createApp(App).mount('#app')
