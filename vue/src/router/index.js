import { createRouter, createWebHistory } from 'vue-router'

const routes = [
  {
    path: '/',
    name: 'overview',
    component: () => import('@swanlab-vue/views/home/HomeView.vue')
  },
  {
    path: '/experiment/:experimentId',
    name: 'experiment',
    component: () => import('@swanlab-vue/views/experiment/ExperimentView.vue')
  },
  {
    path: '/help',
    name: 'help',
    component: () => import('@swanlab-vue/views/help/HelpView.vue')
  }
]

const router = createRouter({
  history: createWebHistory(),
  base: '/',
  routes
})

export default router
