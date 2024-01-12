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
 * 以合理的格式展示数据
 *
 * 常规数据区间为：[1e-7, 1e13)
 * 1. 区间内，返回非科学计数法格式的数据
 * 2. 超出区间，统一使用科学计数法表示
 * 3. 科学计数法底数最多保留 digits 位小数
 * 4. 非科学计数法格式的数据，小数点后，保留到第一个有效数字及其后三位
 * 5. 可选参数 pure 若设置为 true，则在上面四条的基础上，删除无意义的0
 *    - 非科学计数法时，小数点后多余的 0
 *    - 科学计数法时，保留位后多余的 0
 *
 * @param {*} value
 * @param {*} pure
 * @returns
 */
const max = 1e13
const min = 1e-7
const digits = 4
export const formatNumber = (value, pure) => {
  // 检测负数
  const negative = value >= 0 ? false : true
  if (negative) value = Math.abs(value)

  // 如果在区间内，使用普通计数法，只需要保留一定精度即可
  if (value >= min && value < max) {
    value = value.toFixed(13)

    if (pure) {
      value = value.toFixed(digits).replace(/\.?0+$/, '')
    }
  } else {
    value = value.toExponential(digits)
  }

  // 最后处理和返回的时候需要保证是字符串格式，不然，如果符合js自动转化的条件，会自动变成科学计数法
  if (negative) value = '-' + value
  return value.toString()
}
