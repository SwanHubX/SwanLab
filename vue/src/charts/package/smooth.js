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

// 高斯函数
const gaussian = (x, sigma) => {
  return Math.exp((-x * x) / (2 * sigma * sigma))
}

// 高斯平均
const gaussianAverage = (data, param) => {}

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
  if (method.id !== 0) {
    return true
  }
}
