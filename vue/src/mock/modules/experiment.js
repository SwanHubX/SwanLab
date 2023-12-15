import { uuid } from '@swanlab-vue/utils/common'
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
    lossData.push({
      create_date: '1213421312',
      data: lossValue
    })
  }

  return lossData
}

/**
 * 生成随机的RGB值
 * @returns "#RRGGBB" 的形式的字符串
 */
function getRandomHexColor() {
  // 生成随机的 RGB 值
  var r = Math.floor(Math.random() * 256)
  var g = Math.floor(Math.random() * 256)
  var b = Math.floor(Math.random() * 256)

  // 将 RGB 转换为十六进制，并格式化为 "#RRGGBB" 的形式
  var hexColor = '#' + componentToHex(r) + componentToHex(g) + componentToHex(b)

  return hexColor
}

function componentToHex(c) {
  var hex = c.toString(16)
  return hex.length == 1 ? '0' + hex : hex
}

export default [
  // 拿到指定实验的信息
  {
    url: prefix + '/:experiment_id',
    method: 'get',
    response: (experiment_id) => {
      return {
        code: 0,
        message: 'success',
        data: {
          experiment_id: experiment_id,
          name: uuid() + '-' + experiment_id,
          status: 1,
          description: 'this is a test experiment',
          config: {
            learning_rate: 0.01,
            epochs: 10000
          },
          system: {
            hostname: 'caizirui-Pro.local',
            os: 'macOS-14.2-arm64-arm-64bit',
            python: '3.11.5'
          },
          argv: ['E:\\BlackSwan\\swanlab\\test\\test_database_create.py'],
          index: 0,
          tags: ['loss'],
          color: getRandomHexColor(),
          create_time: '2023-12-04T04:44:33.026550',
          update_time: '2023-12-04T04:45:02.036137'
        }
      }
    }
  },
  {
    // 拿到指定实验的指定tag的信息
    url: prefix + '/:experiment_id/tag/:tag', //请求地址
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
  },
  {
    // 获取实验摘要总结信息
    url: prefix + '/:experiment_id/summary',
    method: 'get',
    response: () => {
      return {
        code: 0,
        message: 'success',
        data: {
          summaries: [
            ['loss', 1.23123123],
            ['hello', 332.321312],
            ['world', 3.232323]
          ]
        }
      }
    }
  }
]
