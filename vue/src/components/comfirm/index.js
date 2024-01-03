import SLComfirmVue from './SLComfirm.vue'

const confirmObj = {
  // 用于存储组件实例
  _r: null,
  get _ref() {
    if (!confirmObj._r) console.error('Message组件未注册')
    return confirmObj._r.value
  }
}

/**
 * 弹出确认框
 *
 * @param {string} title - 确认框的标题，如果不传默认为"Are you sure?"。
 * @param {string} content - 确认框的内容, 如果不传默认为"Are you sure you want to do this?"。
 * @param {Object} [config] - 模态对话框的配置选项
 * @param {Object} [config.buttonText] - 按钮的文本, 分为confirm和cancel两个按钮, 如果不传默认为"Yes, I confirm"和"Cancel"。
 * @param {string} [config.buttonText.confirm] - 确认按钮的文本, 如果不传默认为"Yes, I confirm"。
 * @param {string} [config.buttonText.cancel] - 取消按钮的文本, 如果不传默认为"Cancel"。
 *
 * @returns {Promise} - 返回一个Promise对象，当用户点击确认按钮时，Promise会resolve，否则会reject。
 */
export const confirm = (text, content, config = {}) => {
  return new Promise((resolve, reject) => {
    confirmObj._ref.show(text, content, config, resolve, reject)
  })
}

export const installConfirm = (ref) => {
  console.log(ref)
  confirmObj._r = ref
}

export const SLComfirm = SLComfirmVue
