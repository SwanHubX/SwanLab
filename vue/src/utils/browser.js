/** 保存一些和浏览器相关的东西 */
// 1. 判断是否是苹果公司的设备
export const isApple = /iphone|ipad|ipod|mac/i.test(navigator.userAgent)
// 2. 当前设备是否在移动端
export const isMobile = /Mobile/.test(navigator.userAgent)

/**
 * 设置浏览器不可点击，用于异步请求时
 * @param { Function } asyncFunction 异步函数，传入一个请求函数，当请求函数执行完毕后，会自动恢复浏览器的点击事件
 * @returns
 */
export const setBrowserUnClickable = async (asyncFunction) => {
  document.body.style.pointerEvents = 'none'
  return await asyncFunction().finally(() => {
    document.body.style.pointerEvents = 'auto'
  })
}

/**
 * 复制文本到剪贴板，异步函数，使用现代浏览器的API
 * @param { String } text 要复制的文本
 * @param { Function } callback 回调函数 如果执行，回传复制成功与否的状态
 */
export const copyTextToClipboard = async (text, callback) => {
  /**
   * 使用旧的API复制文本到剪贴板，兼容性更好
   * @param { string } text 要复制的文本
   * @param { Function } callback 回调函数 如果执行，回传复制成功与否的状态
   */
  function fallbackCopyToClipboard(text, callback) {
    try {
      const textArea = document.createElement('textarea')
      textArea.value = text
      document.body.appendChild(textArea)
      textArea.select()
      document.execCommand('copy')
      document.body.removeChild(textArea)
      callback(true)
    } catch (e) {
      callback(false)
    }
  }
  let status = true
  try {
    await navigator.clipboard.writeText(text)
  } catch (e) {
    status = false
  }

  if (!status) {
    fallbackCopyToClipboard(text, (success) => {
      status = success
    })
  }
  callback && callback(status)
}

/**
 * 在浏览器主线程添加一个任务
 */
export const addTaskToBrowserMainThread = (task) => {
  setTimeout(task, 0)
}

/**
 * 生成一个uuid，用于标识一个唯一的id
 * @param { Number } now 传入一个时间戳，用于生成uuid，如果不传入，默认为当前时间戳
 * @returns { String } 生成的uuid
 */
export const uuid = (now = Math.round(new Date() / 1000)) => {
  return now.toString(36) + Math.random().toString(36).slice(-4)
}
