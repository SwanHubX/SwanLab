import { defineStore } from 'pinia'
import { ref, computed } from 'vue'

// 这是一个全局的 store，用于存储当前项目的信息，以及一些前端需要的信息
export const useProjectStore = defineStore('project', () => {
  /** state */
  // 当前项目信息
  const project = ref()
  // 当前主题暗还是亮
  const dark = ref(false)

  /** getter */
  const name = computed(() => (project.value?.name ? project.value.name : 'my-machine-learning-project'))
  const description = computed(() => (project.value?.description ? project.value.description : ''))
  const experiments = computed(() => project.value?.experiments)
  const sum = computed(() => project.value?.experiments.length)
  const createTime = computed(() => project.value?.create_time)
  const updateTime = computed(() => project.value?.update_time)
  const logdir = computed(() => project.value?.logdir)
  // 色盘
  const colors = computed(() => {
    const key = dark.value ? 'dark' : 'light'
    return project.value?.colors[key]
  })
  // 当前主题下的所有实验颜色映射关系
  const colorMap = computed(() => {
    const key = dark.value ? 'dark' : 'light'
    const colors = {}
    project.value?.experiments.forEach((e) => {
      colors[e.name] = e[key]
    })
    return colors
  })
  // 当前所有实验可见映射关系
  const showMap = computed(() => {
    const showMap = {}
    project.value?.experiments.forEach((e) => {
      showMap[e.name] = e.show
    })
    return showMap
  })

  /** action */
  const setProject = (p) => {
    p.experiments.reverse()
    project.value = p
  }
  /**
   * 通过实验名获取实验 run_id
   */
  const getExpRunIdByName = (name) => {
    const experiment = experiments.value.find((e) => e.name === name)
    return experiment?.run_id
  }
  /**
   * 删除对应的实验
   * @param {int} id 实验id
   */
  const deleteExperiment = (id) => {
    project.value.experiments = project.value.experiments.filter((e) => e.experiment_id !== id)
  }
  /**
   * 修改项目信息
   * @param {string} name 项目名称
   * @param {string} description 项目描述
   */
  const updateInfo = ({ name, description }) => {
    project.value.name = name
    project.value.description = description
  }
  /**
   * 清空 project
   */
  const clearProject = () => {
    project.value = null
  }
  /**
   * 设置指定实验的状态
   * @param { number } id 实验id
   * @param { string } status 实验状态
   * @param { number } finish_time 实验更新时间
   */
  const setExperimentStatus = (id, status, finish_time) => {
    // 不需要绝对等于
    const experiment = experiments.value.find((e) => e.experiment_id == id)
    experiment.status = status
    if (finish_time) experiment.finish_time = finish_time
  }
  /**
   * 重置实验列表中某一个实验的信息
   * @param {number} id 实验唯一id
   * @param {object} newInfo 涵盖字段为experiment中字段的子集，不需要修改的不必传
   * @param {boolean} overwirte 设置true后，newInfo直接覆盖experiment，请保证newInfo的完整性
   * @returns
   */
  const setExperimentInfo = (id, newInfo, overwirte) => {
    let experiment = experiments.value.find((e) => e.experiment_id == id)
    if (overwirte) {
      experiment = newInfo
      return
    }
    // 遍历 newInfo 中的 key
    for (const key in newInfo) {
      // 如果experiment中没有这个key，说明key不合法
      if (!experiment[key]) continue
      // key 合法，更新该字段
      experiment[key] = newInfo[key]
    }
  }

  /**
   * 修改实验的显示状态，是否显示，二元状态，如果show为1，则改为0，反之亦然
   */
  let changeShowCallback = null
  const changeExperimentShow = (id) => {
    const experiment = experiments.value.find((e) => e.experiment_id == id)
    experiment.show = experiment.show ? 0 : 1
    changeShowCallback && changeShowCallback(id, experiment.show)
    return experiment.show
  }

  const registerChangeShowCallback = (callback) => {
    changeShowCallback = callback
  }
  const destoryChangeShowCallback = () => {
    changeShowCallback = null
  }

  return {
    sum,
    name,
    description,
    experiments,
    createTime,
    updateTime,
    logdir,
    colors,
    colorMap,
    showMap,
    setProject,
    getExpRunIdByName,
    deleteExperiment,
    updateInfo,
    clearProject,
    setExperimentStatus,
    setExperimentInfo,
    changeExperimentShow,
    registerChangeShowCallback,
    destoryChangeShowCallback
  }
})
