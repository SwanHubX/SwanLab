const prefix = '/api/v1/project'

export default [
  // 用户登录接口
  {
    url: prefix + '/experiments', //请求地址
    method: 'get', //请求方式
    response: () => {
      return {
        code: 0,
        message: 'success',
        data: {
          _sum: 2,
          experiments: [
            {
              experiment_id: 1,
              name: 'verdant-maple-1',
              status: 1,
              description: 'this is a test experiment',
              config: {
                learning_rate: 0.01,
                epochs: 10000
              },
              argv: ['E:\\BlackSwan\\swanlab\\test\\test_database_create.py'],
              index: 1,
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
              create_time: '2023-12-04T04:44:33.026550',
              update_time: '2023-12-04T04:45:02.036137'
            },
            {
              experiment_id: 2,
              name: 'abundant-pine-2',
              status: 1,
              description: 'this is a test experiment',
              config: {
                learning_rate: 0.01,
                epochs: 10000
              },
              argv: ['E:\\BlackSwan\\swanlab\\test\\test_database_create.py'],
              index: 2,
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
              create_time: '2023-12-04T12:19:24.329606',
              update_time: '2023-12-04T12:19:51.648795'
            }
          ],
          create_time: '2023-12-04T04:44:33.026550',
          update_time: '2023-12-04T04:44:33.026550'
        }
      }
    }
  }
]
