export default [
  // 用户登录接口
  {
    url: '/mock/api/test', //请求地址
    method: 'get', //请求方式
    response: () => {
      return {
        code: 200,
        msg: 'ok',
        data: Math.floor(Math.random() * 30) + 1
      }
    },
  },
]
