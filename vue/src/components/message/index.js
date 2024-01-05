// ---------------------------------- Message Component,在此处写入组件的使用、注册、导出等信息 ----------------------------------

import SLMessagesVue from './SLMessages.vue'
export const message = {
  // 用于存储组件实例
  _r: null,
  get _ref() {
    if (!message._r) console.error('Message组件未注册')
    return message._r.value
  },
  /**
   * @description 显示错误消息
   * @param {String} text 消息文本
   * @param {Number} delay 延迟关闭时间
   * @param {Function} onClose 关闭时的回调函数
   */
  error: (text, delay, onClose) => {
    message._ref.add(text, 'error', delay, onClose)
  },
  /**
   * @description 显示成功消息
   * @param {String} text 消息文本
   * @param {Number} delay 延迟关闭时间
   * @param {Function} onClose 关闭时的回调函数
   */
  success: (text, delay, onClose) => {
    message._ref.add(text, 'success', delay, onClose)
  },
  /**
   * @description 显示警告消息
   * @param {String} text 消息文本
   * @param {Number} delay 延迟关闭时间
   * @param {Function} onClose 关闭时的回调函数
   */
  warning: (text, delay, onClose) => {
    message._ref.add(text, 'warning', delay, onClose)
  },
  /**
   * 清空所有消息
   */
  clear: () => {
    message._ref.clear()
  }
}

export const installMessage = (ref) => {
  message._r = ref
}

export const SLMessages = SLMessagesVue
