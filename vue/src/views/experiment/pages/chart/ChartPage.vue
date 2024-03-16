<template>
  <!-- 图表容器 -->
  <div class="chart-page" v-if="status === 'success'">
    <ChartsDashboard
      :groups="groups"
      :update-chart-status="updateChartStatus"
      :update-namespace-status="updateNamespaceStatus"
      :media="media"
      :default-color="defaultColor"
      :get-color="getColor"
      :subscribe="on"
      :unsubscribe="off"
    />
    <!-- 图表不存在 -->
    <p class="font-semibold pt-5 text-center" v-if="groups.length === 0">Empty Chart</p>
  </div>
</template>

<script setup>
/**
 * @description: 实验图表页，在本处将完成图表布局的初始化，注入订阅依赖等内容
 * @file: ChartPage.vue
 * @since: 2023-12-25 15:34:51
 **/
import { useExperimentStore, useProjectStore } from '@swanlab-vue/store'
import http from '@swanlab-vue/api/http'
import { ref } from 'vue'
import { useRoute } from 'vue-router'
import { onUnmounted } from 'vue'
import { updateChartStatus, updateNamespaceStatus, media } from '@swanlab-vue/api/chart'
import ChartsDashboard from '@swanlab-vue/charts/ChartsDashboard.vue'
const projectStore = useProjectStore()
const experimentStore = useExperimentStore()
const route = useRoute()

// ---------------------------------- 颜色配置，注入色盘 ----------------------------------
const defaultColor = projectStore.colors[0]
const getColor = (() => {
  // 遍历所有实验，实验名称为key，实验颜色为value
  return (_, index) => {
    return projectStore.colors[index]
  }
})()

// ---------------------------------- 主函数：获取图表配置信息，渲染图表布局 ----------------------------------
// 基于返回的namespcaes和charts，生成一个映射关系,称之为cnMap
// key是chart.id, value是chart.id在charts中的索引
let cnMap = new Map()
let charts = []
const groups = ref([])
const status = ref('initing')
;(async function () {
  const { data } = await http.get(`/experiment/${experimentStore.id}/chart`)
  // 遍历完了，此时cnMap拿到了,模版部分可以依据这些东西渲染了
  parseCharts(data)
  // console.log(cnMap)
})().then(() => {
  // console.log(eventEmitter)
  // status在start函数中被设置为success
  eventEmitter.start()
  // 设置状态为success
  status.value = 'success'
})

// ---------------------------------- 根据charts生成对应配置 ----------------------------------
/**
 * 解析charts，生成groups和cnMap
 */
const parseCharts = (data) => {
  charts = data.charts || []
  groups.value = data.namespaces || []
  // 添加一个临时array，用于减少遍历次数
  const chartsIdArray = charts.map((chart) => {
    // 为每一个chart生成一个cid
    chart._cid = cid(chart.id)
    // 事实上这块需要注意的是，如果chart.error存在，那么chart.source将不会被渲染
    // 但是为了保持代码平衡，这里还是需要将chart.source添加到eventEmitter中
    // 这将在下面第一次请求数据的时候因为错误被停止
    // 往轮询器中添加当前chart的source，即使重复也没事，因为轮询器中是set，默认去重
    eventEmitter.addSource(chart.source)
    console.info('addSource', chart.source, eventEmitter._sources)
    // 返回id
    return chart.id
  })
  // 标志，判断groups是否找到对应的chart
  let found = false
  // 遍历groups
  const temp = groups.value
  groups.value.forEach((group) => {
    found = false
    // 第二次调用以后，id是chart配置，所以id.id
    group.charts.forEach((id) => {
      // 如果id是Object，说明是一个chart配置，取id.id
      if (typeof id === 'object') id = id.id
      // 遍历charts，如果id和charts.charts_id相同，cnMap添加一个key
      for (let i = 0; i < chartsIdArray.length; i++) {
        const chartId = chartsIdArray[i]
        if (chartId === id) {
          found = true
          cnMap.set(chartId, i)
          break
        }
      }
      if (!found) {
        console.error('chart.id: ', id)
        console.error('charts: ', charts)
        console.error('groups: ', groups)
        throw new Error('Charts and groups cannot correspond.')
      }
    })
  })
  // 最终传入的groups，是将charts中的chart_id替换为chart配置
  temp.forEach((group) => {
    group.charts = group.charts.map((id) => {
      if (typeof id === 'object') return charts[cnMap.get(id.id)]
      return charts[cnMap.get(id)]
    })
  })
  groups.value = temp
}

// ---------------------------------- 发布、订阅模型 ----------------------------------
// 请求映射表，反映映射状态
class RequestMap {
  constructor() {
    // 映射表，key:tag ; value:'pending' 'success' 'error'
    this._map = new Map()
  }
  /**
   * 判断是否可以发起请求,只有第一次请求和上一次请求为success的时候能发起请求
   */
  canRequest(key) {
    // 如果key不存在或者key对应的value不为pending，返回true
    // console.log(!this._map.has(key) || this._map.get(key) === 'success')
    // console.log(this._map)
    return !this._map.has(key) || this._map.get(key) === 'success'
  }

  /**
   * 将某个key的状态设置为等待中
   * @param { string } key 对应tag
   */
  setPending(key) {
    this._map.set(key, 'pending')
  }

  /**
   * 将某个key的状态设置为成功
   * @param { string } key 对应tag
   */
  setSuccess(key) {
    this._map.set(key, 'success')
  }

  /**
   * 将某个key的状态设置为失败
   * @param { string } key 对应tag
   */
  setError(key) {
    this._map.set(key, 'error')
  }

  /**
   * 获取某个key的状态
   */
  get(key) {
    return this._map.get(key)
  }
}

// 轮询映射表，提供轮询间隔的设置、依据id清除轮询等功能
class IntervalMap {
  constructor() {
    // 映射表，key:tag ; value:timer
    this._map = new Map()
    // 轮询间隔，key:tag ; value:interval
    this._intervalMap = new Map()
    // 回调函数，key:tag ; value:callback
    this._callbackMap = new Map()
  }
  /**
   * 输入长度，设置轮询间隔，如果轮询间隔没有变化，不会重新设置
   * 1. 数据量小于100，每1秒请求一次
   * 2. 数据量大于100，小于500，每1秒请求一次
   * 3. 数据量大于500，小于1000，每2秒请求一次
   * 4. 之后每增加500个数据，延时增加1秒，最大延时为8秒
   *
   * @param { string } key 对应tag
   * @param { number } length 数据量，依据数据长度设置轮询间隔
   * @param { Function } callback 回调函数，当轮询间隔发生变化时执行
   */
  setInterval(key, length, callback = undefined) {
    // 如果length不存在，直接返回
    if (length === undefined) return
    // 如果callback不存在，且key对应的callback也不存在，直接返回
    if (!callback && !this._callbackMap.get(key)) return
    // console.warn('setInterval', key, length)
    // 计算时间
    let interval = 1000
    if (length < 100) interval = 1000
    else if (length < 500) interval = 1000
    else if (length < 1000) interval = 2000
    else {
      const n = Math.floor((length - 1000) / 500)
      interval = Math.min(8000, 2000 + n * 1000)
    }
    // console.log('setInterval', key, interval)
    if (this._intervalMap.get(key) === interval) return // 如果轮询间隔没有变化，不会重新设置
    // 此时重新设置轮询间隔，并且清除旧的轮询
    this.clear(key)
    // 重新设置轮询
    // console.warn('重新设置轮询', key, interval)
    // 拿到回调函数
    callback = callback || this._callbackMap.get(key)
    this._map.set(key, setInterval(callback, interval))
    this._callbackMap.set(key, callback)
    // 设置轮询间隔
    this._intervalMap.set(key, interval)
  }

  /**
   * 清除轮询
   * @param { string } key 对应tag
   */
  clear(key) {
    const timer = this._map.get(key)
    if (timer) clearInterval(timer)
    this._map.delete(key)
    this._intervalMap.delete(key)
    this._callbackMap.delete(key)
  }

  /**
   * 清除所有轮询
   */
  clearAll() {
    console.log('清除所有轮询')
    this._map.forEach((timer) => {
      console.log('clear', timer)
      clearInterval(timer)
    })
    this._map.clear()
    this._intervalMap.clear()
  }
}

// 事件发射器
class EventEmitter {
  constructor() {
    // 源列表，用于完成数据请求，去重
    this._sources = new Set()
    // 基于源列表形成的源对象，key为源列表的元素，value为请求到的数据
    this._sourceMap = new Map()
    // 订阅缓存，当_sourceMap的状态改变时更新相关内容，它与_sourceMap的key相同
    this._distributeMap = new Map()
    // 实验id
    this._experiment_id = experimentStore.id
    // 轮询id
    this.timer = null
    // 轮询id的映射表，key为tag，value为timer
    this._timerMap = new IntervalMap()
    // 请求间隔映射表，当此表中的数据发生变化时，重新设置轮询函数
    this._intervalMap = new Map()
    // 请求状态映射
    this._request = new RequestMap()
    // 是否执行过start函数
    this._started = false
  }

  /**
   * 更新_sourceMap中的内容，并且完成key的分发
   * @param { string } key 源名称，对应_sourceMap的key
   * @param { Object } data 源数据，会将旧的数据覆盖
   * @param { Object } error 错误信息，如果存在，data将被设置为null
   */
  setSourceData(key, data, error) {
    if (data) {
      this._sourceMap.set(key, data)
    }
    console.log('setSourceData', key, data, error)
    // 遍历对应key的_distributeMap，执行回调函数
    const callbacks = this._distributeMap.get(key)
    // 回调函数必须全部执行完毕后，才能将请求状态设置为success或者error
    if (!callbacks) return this.setRequestStatus(key, error)
    // callbacks存在，则是一个list，每个元素是一个回调函数，回调函数返回一个promise，当所有promise执行完毕后，设置请求状态
    Promise.all(Array.from(callbacks.values()).map((callback) => callback(key, data, error))).then(() => {
      console.log('all callbacks are executed')
      this.setRequestStatus(key, error)
      this._timerMap.setInterval(key, data?.list?.length)
    })
  }

  /**
   * 将source添加到模型的源列表中
   * @param { Array } source 源列表，对应每个chart的source
   */
  addSource(source) {
    source = source || []
    source.map((s) => {
      // 如果还没有执行过start函数，直接添加到源列表中
      if (!this._started) return this._sources.add(s)
      // 判断这个数据是否存在，即使set已经去重，但是需要开启轮询
      if (this._sources.has(s)) return
      // 添加到源列表中
      this._sources.add(s)
      // 发起请求
      console.warn('添加源并开启轮询: ', s)
      this._getTagRequestInterval(s)
    })
  }

  /**
   * 启动事件发射器
   * 执行此函数，前端开始获取tag数据，并且依据实验状态判断是否关闭轮询等
   *
   */
  start() {
    // 遍历源列表，请求数据
    const promises = []
    this._sources.forEach((tag) => {
      // promise all如果出现错误，会直接reject，不会执行后面的，所以这里不用它,使用for of
      promises.push(this._getSoureceData(tag))
    })
    // 等待第一次请求完成，然后设置轮询
    Promise.all(promises).then(() => {
      // 如果实验状态不是0，停止轮询
      if (!experimentStore.isRunning) return console.log('stop, experiment status is not 0')
      // 遍历源列表，请求数据
      this._sources.forEach((tag) => {
        // 设置轮询
        this._getTagRequestInterval(tag)
      })
      this._started = true
    })
  }

  /**
   * 运行此函数代表必然开启轮询
   * @param { string } tag 源名称
   */
  _getTagRequestInterval(tag) {
    // console.log('设置轮询', tag)
    this._timerMap.setInterval(tag, this._sourceMap.get(tag)?.list.length || 0, () => {
      // 在此处取出pinia中的charts配置，重新设置
      experimentStore.charts && parseCharts(experimentStore.charts)
      // 判断当前实验id和路由中的实验id是否相同，如果不同，停止轮询
      if (Number(this._experiment_id) !== Number(route.params.experimentId)) {
        this._timerMap.clearAll()
        return console.log('stop, experiment id is not equal to route id')
      }
      this._getSoureceData(tag)
      // 如果实验状态不是0，停止轮询
      if (!experimentStore.isRunning) {
        this._timerMap.clearAll()
        return console.log('stop, experiment status is not 0')
      }
    })
  }

  /**
   * 根据source获取数据，此时source又被称为tag
   * @param { string } source 源名称
   */
  _getSoureceData(tag) {
    // 如果无法请求，直接返回
    if (!this._request.canRequest(tag)) return console.log('abandon')
    // 在请求之前将_request设置为pending
    this._request.setPending(tag)
    return http
      .get(`/experiment/${this._experiment_id}/tag/${tag}`)
      .then((res) => {
        // 更新_sourceMap
        console.log('更新')
        this.setSourceData(tag, res.data, null)
      })
      .catch((err) => {
        // 更新_sourceMap
        this.setSourceData(tag, null, err)
      })
  }

  /**
   * 设置请求状态
   * @param { string } tag  源名称
   * @param { object } error 错误信息，如果存在，设置为error，否则设置为success
   */
  setRequestStatus(tag, error) {
    if (error) this._request.setError(tag)
    else this._request.setSuccess(tag)
  }

  /**
   * 订阅数据，当数据更新时，执行回调函数
   * @param { array } tags 源列表，对应chart的source，可以是一个或多个
   * @param { string } cid 订阅者id，对应chart的_cid
   * @param { Function } callback 回调函数，当数据更新时执行订阅者的回调函数，这应该是一个Promise
   * 回调函数接收三个参数，分别是tag、data和err
   */
  $on(tags, cid, callback) {
    console.log('订阅', tags, cid)
    tags = tags || []
    tags.map((tag) => {
      // 拿到已经订阅的回调函数
      const callbacks = this._distributeMap.get(tag)
      // 如果已订阅存在，直接添加
      if (callbacks) callbacks.set(cid, callback)
      // 否则，创建一个新的map，添加
      else this._distributeMap.set(tag, new Map([[cid, callback]]))
      // 如果这个tag的请求并不存在，发起请求
      if (this._request.canRequest(tag)) {
        console.warn('数据不存在，发起请求')
        return this._getSoureceData(tag)
      }
      // 如果这个tag的请求已经失败，执行回调函数
      if (this._request.get(tag) === 'error') {
        console.warn('数据请求失败，执行回调函数')
        return callback(tag, null, true)
      }
    })
  }

  /**
   * 取消订阅，当chart被销毁时，需要手动取消订阅
   * @param { string } tags 源列表，对应chart的source，可以是一个或多个
   * @param { string } cid 订阅者id，对应chart的_cid
   */
  $off(tags, cid) {
    if (!this?._distributeMap) return
    tags = tags || []
    tags.map((tag) => {
      // 拿到已经订阅的回调函数
      const callbacks = this._distributeMap.get(tag)
      // 如果已订阅存在，直接添加
      if (callbacks) callbacks.delete(cid)
    })
  }

  /**
   * 销毁订阅模型
   */
  delete() {
    this._timerMap.clearAll()
    // this._sources = null
    // this._sourceMap = null
    // this._distributeMap = null
  }
}

// 初始化订阅模型
const eventEmitter = new EventEmitter()

// 依赖注入，注入订阅函数和取消订阅函数
const on = (tags, cid, callback) => eventEmitter.$on(tags, cid, callback)
const off = (tags, cid) => eventEmitter.$off(tags, cid)

// 当前组件销毁时，取消相关行为
onUnmounted(() => {
  eventEmitter.delete()
})

// ---------------------------------- 工具函数：辅助渲染 ----------------------------------
/**
 * 依据chart_id生成cid，这是一个唯一id，用于eventEmitter中的订阅id，也是dom对象的id
 * @param { string } chart_id
 */
const cid = (chart_id) => {
  return chart_id
}
</script>

<style lang="scss" scoped>
.chart-page {
  @apply bg-higher min-h-[calc(100vh-163px)];
}
</style>
