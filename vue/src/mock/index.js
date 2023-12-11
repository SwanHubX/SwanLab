/** 配置前端 mock */
import { createProdMockServer } from 'vite-plugin-mock/es/createProdMockServer'
import project from './modules/project'
import experiment from './modules/experiment'

export function setupProdMockServer() {
  createProdMockServer([...project, ...experiment])
}
