<template>
  <!-- 标题 -->
  <div class="flex items-center gap-2 py-5 px-8 border-b">
    <div class="w-5 h-5 rounded-full" :style="{ backgroundColor: experimentStore.color }"></div>
    <h1 class="text-lg font-semibold">{{ experimentStore.name }}</h1>
    <SLStatusLabel :status="experimentStore.status" />
  </div>
  <!-- 图表容器 -->
  <template v-if="status === 'success'">
    <ChartsContainer
      v-for="group in groups"
      :key="group.namespace"
      :label="getGroupName(group.namespace)"
      :count="group.charts.length"
    >
      <ChartContainer :chart="charts[cnMap.get(chartID)]" v-for="chartID in group.charts" :key="chartID" />
    </ChartsContainer>
  </template>
  <EmptyExperiment v-once v-else-if="status !== 'initing'" />
</template>

<script setup>
/**
 * @description: 实验图表页，在本处将完成图表布局的初始化，注入订阅依赖等内容
 * @file: ChartPage.vue
 * @since: 2023-12-25 15:34:51
 **/
import SLStatusLabel from '@swanlab-vue/components/SLStatusLabel.vue'
import { useExperimentStroe } from '@swanlab-vue/store'
import http from '@swanlab-vue/api/http'
import { ref, provide } from 'vue'
import ChartsContainer from './components/ChartsContainer.vue'
import ChartContainer from './components/ChartContainer.vue'
import EmptyExperiment from './components/EmptyExperiment.vue'
import { t } from '@swanlab-vue/i18n'
import { useRoute } from 'vue-router'
import { onUnmounted } from 'vue'
const experimentStore = useExperimentStroe()
const route = useRoute()

// ---------------------------------- 主函数：获取图表配置信息，渲染图表布局 ----------------------------------
// 基于返回的namespcaes和charts，生成一个映射关系,称之为cnMap
// key是chart_id, value是chart_id在charts中的索引
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
  status.value = 'success'
  // 在完成以后，开始请求数据
  // console.log(eventEmitter)
  eventEmitter.start()
})

// ---------------------------------- 根据charts生成对应配置 ----------------------------------

const parseCharts = (data) => {
  charts = data.charts || []
  groups.value = data.namespaces || []
  // 添加一个临时array，用于减少遍历次数
  const chartsIdArray = charts.map((chart) => {
    // 为每一个chart生成一个cid
    chart._cid = cid(chart.chart_id)
    eventEmitter.addSource(chart.source)
    // 返回id
    return chart.chart_id
  })
  // 标志，判断groups是否找到对应的chart
  let found = false
  // 遍历groups
  const temp = groups.value
  groups.value.forEach((group) => {
    found = false
    group.charts.forEach((id) => {
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
        console.error('charts: ', charts)
        console.error('groups: ', groups)
        throw new Error('Charts and groups cannot correspond.')
      }
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
    // 请求状态映射
    this._request = new RequestMap()
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
    // 遍历对应key的_distributeMap，执行回调函数
    const callbacks = this._distributeMap.get(key)
    // 回调函数必须全部执行完毕后，才能将请求状态设置为success或者error
    if (!callbacks) return this.setRequestStatus(key, error)
    // callbacks存在，则是一个list，每个元素是一个回调函数，回调函数返回一个promise，当所有promise执行完毕后，设置请求状态
    Promise.all(Array.from(callbacks.values()).map((callback) => callback(key, data, error))).then(() => {
      console.log('all callbacks are executed')
      this.setRequestStatus(key, error)
    })
  }

  /**
   * 将source添加到模型的源列表中
   * @param { Array } source 源列表，对应每个chart的source
   */
  addSource(source) {
    source = source || []
    source.map((s) => this._sources.add(s))
  }

  /**
   * 启动事件发射器
   * 执行此函数，前端开始获取tag数据，并且依据实验状态判断是否关闭轮询等
   *
   */
  start() {
    // 第一次延时1秒执行，后面每隔n秒执行一次
    const n = 1
    // 遍历源列表，请求数据
    this._sources.forEach((tag) => {
      // promise all如果出现错误，会直接reject，不会执行后面的，所以这里不用它,使用for of
      this._getSoureceData(tag)
    })
    if (experimentStore.isRunning)
      this.timer = setInterval(async () => {
        // 在此处取出pinia中的charts配置，重新设置
        experimentStore.charts && parseCharts(experimentStore.charts)
        // 判断当前实验id和路由中的实验id是否相同，如果不同，停止轮询
        if (Number(this._experiment_id) !== Number(route.params.experimentId)) {
          clearInterval(this.timer)
          return console.log('stop, experiment id is not equal to route id')
        }
        // 遍历源列表，请求数据
        this._sources.forEach((tag) => {
          // promise all如果出现错误，会直接reject，不会执行后面的，所以这里不用它,使用for of
          this._getSoureceData(tag)
        })
        // 如果实验状态不是0，停止轮询
        if (experimentStore.status !== 0) {
          clearInterval(this.timer)
          return console.log('stop, experiment status is not 0')
        }
      }, n * 1000)
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
    http
      .get(`/experiment/${this._experiment_id}/tag/${tag}`)
      .then((res) => {
        // 更新_sourceMap
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
    tags = tags || []
    tags.map((tag) => {
      // 拿到已经订阅的回调函数
      const callbacks = this._distributeMap.get(tag)
      // 如果已订阅存在，直接添加
      if (callbacks) callbacks.set(cid, callback)
      // 否则，创建一个新的map，添加
      else this._distributeMap.set(tag, new Map([[cid, callback]]))
      // 如果这个tag的请求并不存在，发起请求
      if (this._request.canRequest(tag)) return this._getSoureceData(tag)
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
    clearInterval(this.timer)
    // this._sources = null
    // this._sourceMap = null
    // this._distributeMap = null
  }
}

// 初始化订阅模型
const eventEmitter = new EventEmitter()

// 依赖注入，注入订阅函数和取消订阅函数
provide('$on', (tags, cid, callback) => eventEmitter.$on(tags, cid, callback))
provide('$off', (tags, cid) => eventEmitter.$off(tags, cid))

// 当前组件销毁时，取消相关行为
onUnmounted(() => {
  eventEmitter.delete()
})

// ---------------------------------- 工具函数：辅助渲染 ----------------------------------

/**
 * 将组名进行一些翻译
 * @param { string } namespace 组名
 */
const getGroupName = (namespace) => {
  console.log(namespace)
  if (namespace === 'default') return t('experiment.chart.label.default')
  else return namespace
}

/**
 * 依据chart_id生成cid，这是一个唯一id，用于eventEmitter中的订阅id，也是dom对象的id
 * @param { string } chart_id
 */
const cid = (chart_id) => {
  return chart_id
}
</script>

<style lang="scss" scoped></style>
