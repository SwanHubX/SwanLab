/**
 * 数据、事件处理等公用函数
 */

/**
 * 基础防抖函数，用于减少函数的调用频率，接收两个参数，第一个是要防抖的函数，第二个是延迟的时间
 * @param { Function } func 函数
 * @param { string } delay 延迟时间，单位毫秒
 * @returns { Function } 返回一个新的函数
 */
export function debounce(func, delay) {
  let timer = null
  return function (...args) {
    if (timer) clearTimeout(timer)
    timer = setTimeout(() => {
      func.apply(this, args)
    }, delay)
  }
}

/**
 * 生成一个uuid，用于标识一个唯一的id
 * @param { Number } now 传入一个时间戳，用于生成uuid，如果不传入，默认为当前时间戳
 * @returns { String } 生成的uuid
 */
export const uuid = (now = Math.round(new Date() / 1000)) => {
  return now.toString(36) + Math.random().toString(36).slice(-4)
}

/**
 * 格式化时间，接收一个时间戳，转化当前浏览器当前时区下的时间，输入时间为UTC时间
 * 例如：上海时区下，2023-12-04T04:44:33.026550转换为"2023/12/04 12:44:33"
 */

export const formatTime = (time) => {
  const date = new Date(time)
  const timezoneOffset = new Date().getTimezoneOffset()
  // console.log(timezoneOffset)
  const localTime = date.getTime() - timezoneOffset * 60 * 1000
  const localDate = new Date(localTime)
  const year = localDate.getFullYear()
  const month = localDate.getMonth() + 1
  const day = localDate.getDate()
  const hour = localDate.getHours()
  const minute = localDate.getMinutes()
  const second = localDate.getSeconds()
  return `${year}/${month}/${day} ${hour}:${minute}:${second}`
}
