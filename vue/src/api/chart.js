/**
 * chart组件相关api，包括获取媒体文件，更改namespace展开状态、chart置顶等
 */
import http from './http'

/**
 * 获取媒体文件，获取blob对象
 * 返回promise，如果成功，返回blob对象，否则返回错误信息
 */
export const media = {
  /**
   * 获取媒体文件，获取blob对象
   * 返回promise，如果成功，返回blob对象，否则返回错误信息
   * @param { string } data 即文件名，即数据中的data字段
   * @param { string } experiment_id 实验id
   * @param { string } tag 数据名称
   * @returns { Promise<Blob> }
   */
  get: (data, experiment_id, tag) => {
    return new Promise((resolve, reject) => {
      http
        .get('/media/' + data, { params: { tag, experiment_id }, responseType: 'blob' })
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
 * 更新命名空间展开状态
 * @param { boolean } opened 是否展开
 * @param { object } namespace 命名空间对象
 * @returns { Promise }
 */
export const updateNamespaceStatus = (opened, namespace) => {
  return http.patch('/namespace/' + namespace.id + '/opened', {
    experiment_id: namespace.experiment_id?.id,
    project_id: namespace.project_id?.id,
    opened
  })
}

/**
 * 更新chart状态，置顶、隐藏或正常显示
 * @param { object } chart 图表对象
 * @param { int } status 状态码，0为正常，1为置顶，-1为隐藏
 * @returns { Promise } Promise对象，最终可返回更新后的图表组织结构
 */
export const updateChartStatus = async (chart, status) => {
  const { data } = await http.patch('/chart/' + chart.id + '/status', {
    status
  })
  return data
}
