import { uuid, formatNumber2SN } from '@swanlab-vue/utils/common'
import { describe, test, expect } from 'vitest'

describe('common.js utils', () => {
  test('generate uuid', () => {
    expect(uuid()).toHaveLength(10)
  })

  /**
   * 对于格式化数据的测试
   * 具体的格式规则可见 /vue/utils/common.js
   * 这里的测试主要针对几个类型：
   * 1. 正负数
   * 2. 区间内的小数(小数点后的精度)
   * 3. 区间内绝对值大于1的数(是否带小数部分)
   * 4. 在第3点的基础上，小数位超出额定范围，抛弃小数
   * 5. 超出区间的极小数
   * 6. 超出区间的极大数
   */

  test('formatNumber2SN test for interval fraction numbers', () => {
    const result1 = formatNumber2SN(0.12345)
    expect(result1).toBe('0.1234')

    const result2 = formatNumber2SN(-0.123456)
    expect(result2).toBe('-0.1234')

    const result3 = formatNumber2SN(-0.000000123456)
    // 注意！这里有四舍五入
    expect(result3).toBe('-0.0000001235')
  })

  test('formatNumber test1 for interval medium sized numbers', () => {
    // 正常的正负数，直接输出
    const result1 = formatNumber2SN(123312312)
    expect(result1).toBe('123312312')

    const result2 = formatNumber2SN(-21312312)
    expect(result2).toBe('-21312312')

    // 对小数部分，只截取到第一个非零数及其后三位，且自动去除无意义的0
    const result3 = formatNumber2SN(1233.0403000001)
    expect(result3).toBe('1233.0403')

    const result4 = formatNumber2SN(-1233.0000123456)
    expect(result4).toBe('-1233.00001234')

    // 小数部分全为0，则只取整数部分
    const result5 = formatNumber2SN(123321.0)
    expect(result5).toBe('123321')
  })

  test('formatNumber test2 for interval medium sized numbers', () => {
    const result1 = formatNumber2SN(12.00000000000001)
    expect(result1).toBe('12')

    const result2 = formatNumber2SN(-12.00000000000001)
    expect(result2).toBe('-12')
  })

  test('formatNumber test for extreme decimals out of range', () => {
    const result1 = formatNumber2SN(0.0000000000000000001234512)
    expect(result1).toBe('1.2345e-19')

    const result2 = formatNumber2SN(-0.0000000000000000001234512)
    expect(result2).toBe('-1.2345e-19')
  })

  test('formatNumber test for maximum number exceeding the interval', () => {
    const result1 = formatNumber2SN(1010000000000123)
    expect(result1).toBe('1.01e+15')

    const result2 = formatNumber2SN(-1010000000000123)
    expect(result2).toBe('-1.01e+15')
  })
})
