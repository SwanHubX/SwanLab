import axios from 'axios'

// axios对象实例
const http = axios.create({
  // baseURL不设置，会自动使用当前域名
  // baseURL: '',
  timeout: 10000
})

http.interceptors.request.use(
  async (req) => {
    // TODO 判断token过期时间
    console.log('[request] ', req.method, req.url, req.data || req.params || '')
    return req
  },
  (error) => {
    console.log('[request error] ', error)
    return Promise.reject(error)
  }
)

// 响应拦截器
http.interceptors.response.use(
  (resp) => {
    // 打印响应信息
    console.log('[response] ', resp.config.url, resp.data)
    return resp.data
  },
  (error) => {
    // 判断连接是否超时
    if (error.code === 'ECONNABORTED' && error.message.indexOf('timeout') !== -1) {
      console.log('[response error] ', 'timeout')
      return Promise.reject(error)
    }
    console.log('[response error] ', error.response?.status, error.response?.data, error)
    return Promise.reject(error)
  }
)

export default http
