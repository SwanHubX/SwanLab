import { createProdMockServer } from 'vite-plugin-mock/es/createProdMockServer'
import MockMethod from './api'

export function setupProdMockServer() {
  createProdMockServer([...MockMethod])
}
