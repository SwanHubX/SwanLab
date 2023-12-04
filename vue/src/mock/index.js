/** 配置前端 mock */
import { createProdMockServer } from 'vite-plugin-mock/es/createProdMockServer'
import test_api from './test'
import project_api from './project'

export function setupProdMockServer() {
  createProdMockServer([
    ...test_api,
    ...project_api
  ])
}
