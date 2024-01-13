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
export const formatNumber = (value, pure = true) => {
  // 检测负数
  const negative = value >= 0 ? false : true
  // 如果是负数，暂时取绝对值运算
  if (negative) value = Math.abs(value)

  // 如果在区间内，使用普通计数法，只需要保留一定精度即可
  if (value >= min && value < max) {
    /**
     * 这里使用 toFixed 的目的是，js 自动转化科学计数法的最小阈值是 e-6，而我们需要 e-7 才转化
     * 为了将 e-7 ~ e-6 之间的数从科学计数法转为十进制，需要转化一下
     * 要求保留小数点后第一个有效数字和其后三位，那么最多需要保留10位精度，例：0.000000123424
     */
    if (value.toString().includes('e')) {
      // 获取最小的精度位
      const min_str = min.toString()
      const min_digits = Number(min_str.substring(min_str.indexOf('-') + 1))
      value = value.toFixed(min_digits + digits - 1)
    }

    // 到这里，可以保证 value 是一个十进制表示的，位于区间内的值，且为 string
    value = value.toString()

    // 找到小数点后第一位有效数字的索引
    const match = value.match(/\.\d*?([1-9])/)
    if (match) {
      // 如果找到匹配，获取 value 有效部分的长度
      // 小数点后第一个非零数在小数点后的索引（match[0]中包含了小数点 => 长度-1）
      const index = match[0].length - 1
      /**
       * 计算符合展示规则的长度
       * 整数部分长度：value.substring(0, value.indexOf('.')).length
       * 小数部分长度：index + digits
       */
      const length = value.substring(0, value.indexOf('.')).length + index + digits
      value = value.substring(0, length)
    } // 没有匹配上说明没有小数部分，直接略过

    /**
     * 默认删除小数点后，无意义的 0
     * 例如：1.01300 => 1.013
     */
    if (pure) {
      if (value.indexOf('.') !== -1) {
        value = value.replace(/(\.\d*?[1-9])?0+$/, '$1')
      }
    }
  }
  // 如果超出范围，统一使用科学计数法
  else {
    value = value.toExponential(digits) // 注意 toExponential 的返回值是 string
    if (pure) value = value.replace(/(\.\d*?[1-9])0+e/, '$1e') // 清除无意义的 0
  }

  // 最后处理和返回的时候需要保证是字符串格式，不然，如果符合js自动转化的条件，会自动变成科学计数法
  if (negative) value = '-' + value
  return value
}
