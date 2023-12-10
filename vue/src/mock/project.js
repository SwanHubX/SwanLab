import { uuid } from '@swanlab-vue/utils/common'
const prefix = '/api/v1/project'

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

const generateExperiments = (num) => {
  let experiments = []
  for (let i = 0; i < num; i++) {
    let status = Math.floor(Math.random() * 3)
    if (status === 2) status = -1
    // 随机生成一个颜色
    experiments.push({
      experiment_id: i,
      name: uuid() + '-' + i,
      status: status,
      description: 'this is a test experiment',
      config: {
        learning_rate: 0.01,
        epochs: 10000
      },
      argv: ['E:\\BlackSwan\\swanlab\\test\\test_database_create.py'],
      index: i,
      tags: [
        {
          tag: 'loss',
          num: 9998
        },
        {
          tag: 'accuracy',
          num: 9998
        }
      ],
      color: getRandomHexColor(),
      create_time: '2023-12-04T04:44:33.026550',
      update_time: '2023-12-04T04:45:02.036137'
    })
  }
  return experiments
}

export default [
  {
    // 列出当前项目下所有实验(实验信息)
    url: prefix, //请求地址
    method: 'get', //请求方式
    response: () => {
      const sum = 20
      const experiments = generateExperiments(sum)
      // console.log('experiments', experiments)
      return {
        code: 0,
        message: 'success',
        data: {
          _sum: sum,
          experiments,
          create_time: '2023-12-04T04:44:33.026550',
          update_time: '2023-12-04T04:44:33.026550'
        }
      }
    }
  }
]
