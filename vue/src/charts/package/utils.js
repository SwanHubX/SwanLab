/**
 * 在此处封装一些图表可共用的函数逻辑
 */
import http from '@swanlab-vue/api/http'
export const media = {
  /**
   * 获取媒体文件，获取blob对象
   * @param { string } data 即文件名，即数据中的data字段
   * @param {*} runId 用于区分不同的实验，即run_id字段
   * @param {*} tag 数据名称
   * @returns { Promise<Blob> }
   */
  get: (data, runId, tag) => {
    return http.get('/media/' + data, { params: { run_id: runId, tag }, responseType: 'blob' })
  }
}

/**
 * 将refrence字段映射为x轴字段和显示的x轴名称
 */
export const refrence2XField = {
  step: { xField: 'index', xTitle: 'Step' }
}
