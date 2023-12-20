import zhCN from './zh-CN'
import en from './en-US'
import { createI18n } from 'vue-i18n'

export const i18n = createI18n({
  locale: 'en', // 设置地区
  legacy: false, // 如果要支持compositionAPI，此项必须设置为false
  globalInjection: true, // 全局注册$t方法
  messages: {
    'zh-CN': zhCN,
    en
  }
})

export const { t, te, tm } = i18n.global
