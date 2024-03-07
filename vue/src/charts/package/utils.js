/**
 * 在此处封装一些图表可共用的函数逻辑
 */
import http from '@swanlab-vue/api/http'
export const media = {
  /**
   * 获取媒体文件，获取blob对象
   * 返回promise，如果成功，返回blob对象，否则返回错误信息
   * @param { string } data 即文件名，即数据中的data字段
   * @param {*} run_id 用于区分不同的实验，即run_id字段
   * @param {*} tag 数据名称
   * @returns { Promise<Blob> }
   */
  get: (data, run_id, tag) => {
    return new Promise((resolve, reject) => {
      http
        .get('/media/' + data, { params: { tag, run_id }, responseType: 'blob' })
        .then((res) => {
          resolve(res)
        })
        .catch((err) => {
          // 如果blob的type为application/json，说明是后端返回的错误信息，需要解析
          if (err.data.type === 'application/json') {
            const reader = new FileReader()
            reader.readAsText(err.data)
            reader.onload = () => {
              reject(JSON.parse(reader.result))
            }
          } else reject(err)
        })
    })
  }
}

/**
 * 将refrence字段映射为x轴字段和显示的x轴名称
 */
export const refrence2XField = {
  step: { xField: 'index', xTitle: 'Step' }
}

export const transparentColor = (color, opacity = 0.1) => {
  // 将十六进制增加透明度通道
  opacity = Math.floor(opacity * 255)
  // opacity转为十六进制
  const opacityHex = opacity.toString(16)
  return color + opacityHex
}
