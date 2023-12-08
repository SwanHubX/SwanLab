/** 配置前端 mock */
import { createProdMockServer } from 'vite-plugin-mock/es/createProdMockServer'
import project from './project'

export function setupProdMockServer() {
  createProdMockServer([...project])
}
