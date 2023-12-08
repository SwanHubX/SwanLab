const prefix = '/api/v1/experiment'

/**
 * 生成逐渐下降的模拟损失数据集，最终趋近于实际损失。
 *
 * @param {number} numSamples - 模拟数据集中的样本数量。
 * @param {number} initialLoss - 初始损失值，数据集的第一个样本的损失。
 * @param {number} finalLoss - 最终损失值，数据集的最后一个样本的损失。
 * @returns {number[]} 包含逐渐下降损失值的数组。
 */
function generateDecreasingLossData(numSamples, initialLoss, finalLoss) {
  const lossData = []

  for (let i = 0; i < numSamples; i++) {
    // 使用指数函数生成逐渐下降的损失值
    const progress = i / (numSamples - 1) // 范围在 [0, 1] 之间
    const lossValue = initialLoss * Math.pow(finalLoss / initialLoss, progress)

    // 将损失值添加到数据集中
    lossData.push(lossValue)
  }

  return lossData
}

export default [
  {
    // 列出当前项目下所有实验(实验信息)
    url: prefix + '/:experiment_id/:tag', //请求地址
    method: 'get', //请求方式
    response: () => {
      return {
        code: 0,
        message: 'success',
        data: {
          list: generateDecreasingLossData(1000, 1.0, 0.1)
        }
      }
    }
  }
]
