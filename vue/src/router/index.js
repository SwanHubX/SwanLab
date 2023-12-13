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
    component: () => import('@swanlab-vue/views/experiment/ExperimentView.vue'),
    redirect: { name: 'experiment_chart' },
    children: [
      {
        path: 'index',
        name: 'experiment_index',
        component: () => import('@swanlab-vue/views/experiment/pages/index/IndexPage.vue')
      },
      {
        path: 'chart',
        name: 'experiment_chart',
        component: () => import('@swanlab-vue/views/experiment/pages/chart/ChartPage.vue')
      },
      {
        path: 'log',
        name: 'experiment_log',
        component: () => import('@swanlab-vue/views/experiment/pages/LogPage.vue')
      }
    ]
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
