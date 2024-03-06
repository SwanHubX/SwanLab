/**
 * 数据平滑处理
 */

// 时间加权均衡
const timeEMA = (data, param) => {
  console.log('tema', data, param)
  // 基础smoothingWeight
  const smoothingWeight = Math.min(Math.sqrt(param || 0), 0.999)
  let debiasWeight = 0
  // const rangeOfX = data[data.length - 1].index - data[0].index
  // 上一个点的y值
  let lastY = data.length > 0 ? 0 : NaN
  return data.map((item, index) => {
    const prevX = index > 0 ? index - 1 : 0
    // VIEWPORT_SCALE scales the result to the chart's x-axis range
    const changeInX = data[index].index - data[prevX].index
    const smoothingWeightAdj = Math.pow(smoothingWeight, changeInX)
    lastY = lastY * smoothingWeightAdj + data[index].data
    debiasWeight = debiasWeight * smoothingWeightAdj + 1
    return { ...item, data: lastY / debiasWeight }
  })
}

// 流动平均
const runningAverage = (data, param) => {
  // 理论上的窗口大小
  const w = Math.floor(param / 2)
  let floorX, ceilX, nowIndex
  let sum = 0
  return data.map((item, index) => {
    floorX = index - w
    ceilX = index + w
    console.log('runningAverage', floorX, ceilX)
    sum = 0
    for (let i = floorX; i <= ceilX; i++) {
      nowIndex = Math.min(Math.max(i, 0), data.length - 1)
      sum += data[nowIndex].data
    }
    return { ...item, data: sum / (ceilX - floorX + 1) }
  })
}

// 创建高斯函数
const createGaussian = () => {
  // 根号2π
  const sqrt2PI = Math.sqrt(2 * Math.PI)
  /**
   * 高斯函数
   * @param { Number } x 自变量
   * @param { Number } stdDev 标准差
   */
  return (x, stdDev) => {
    return Math.exp(-0.5 * Math.pow(x / stdDev, 2)) / (stdDev * sqrt2PI)
  }
}

// 高斯平均
const gaussianAverage = (data, param) => {
  const gaussian = createGaussian()
  // 在代码中正负无穷似乎在大数据上影响算法速度，但是高斯函数在-3到3区间上进行积分似乎已经达到0.99以上，考虑使用此区间即可
  const range = 6
  // 对于每一个值，计算高斯函数的值，只采用-3到3的区间
  let sum = 0
  let weightSum = 0
  return data.map((item, index) => {
    sum = 0
    weightSum = 0
    for (let i = -range; i <= range; i++) {
      if (index + i >= 0 && index + i < data.length) {
        sum += data[index + i].data * gaussian(i, param)
        weightSum += gaussian(i, param)
      }
    }
    return { ...item, data: sum / weightSum }
  })
}

// 平滑算法映射关系
const smoothArithmetics = new Map([
  [0, undefined],
  [1, timeEMA],
  [2, runningAverage],
  [3, gaussianAverage]
])

/**
 * 数据平滑处理，传入数据和平滑方法
 * @param { Array } data 数据
 * @param { Object } method 平滑方法
 * @param { Number } method.id 平滑方法id
 * @param { Number} method.value 平滑方法参数
 * @returns { Array } 平滑后的数据
 */
export default function smooth(data, method) {
  // console.log('smooth', data, method)
  return smoothArithmetics.get(method.id)(data, method.value)
}

/**
 * 判断是否需要平滑，传入平滑方法
 * @param { Object } method 平滑方法
 * @param { Number } method.id 平滑方法id
 * @param { Number} method.value 平滑方法参数
 * @returns { Boolean } 是否需要平滑，true为需要，false为不需要
 */
export const needSmooth = (method) => {
  if (!method.value) {
    return false
  }
  if (method.id !== 0) {
    return true
  }
}
