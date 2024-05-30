<template>
  <div class="flex flex-col justify-between grow" ref="lineChartRef">
    <!-- 图表标题 -->
    <p class="text-center font-semibold select-none">{{ title }}</p>
    <div class="flex flex-col justify-center grow text-dimmer gap-2" v-if="error">
      <SLIcon class="mx-auto h-5 w-5" icon="error" />
      <p class="text-center text-sm">
        {{ $t('chart.charts.line.error', { type: error['data_class'], tag: source[0] }) }}
      </p>
    </div>
    <template v-else>
      <!-- 图表主体 -->
      <LineChartLegend
        :items="legend"
        :key="legend.length"
        @hoverin="(item) => legendHoverin(item, false)"
        @hoverout="(item) => legendHoverout(item, false)"
        v-if="legend && multi"
      />
      <div class="relative">
        <!-- x轴坐标单位 -->
        <p class="absolute right-2 bottom-6 text-xs text-dimmer scale-90 select-none">{{ xTitle }}</p>
        <LineChartTooltip ref="tooltipRef" />
        <div class="overflow-hidden" ref="g2Ref"></div>
      </div>
      <!-- 放大效果 -->
      <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom" esc-exit>
        <p class="text-center mt-4 mb-10 text-2xl font-semibold">{{ title }}</p>
        <LineChartLegend
          :items="legend"
          @hoverin="(item) => legendHoverin(item, true)"
          @hoverout="(item) => legendHoverout(item, true)"
          v-if="legend && multi"
        />
        <div class="relative" ref="g2ZoomRef">
          <p class="absolute right-2 bottom-6 text-xs text-dimmer scale-90 select-none">{{ xTitle }}</p>
          <LineChartTooltip ref="tooltipZoomRef" />
        </div>
      </SLModal>
    </template>
  </div>
</template>

<script setup>
/**
 * @description: 折线图表
 * @file: LineChart.vue
 * @since: 2023-12-25 20:17:19
 **/
import SLModal from '@swanlab-vue/components/SLModal.vue'
import SLIcon from '@swanlab-vue/components/SLIcon.vue'
import { Line, G2 } from '@antv/g2plot'
import * as UTILS from './utils'
import { ref, inject, computed, onUnmounted, provide } from 'vue'
import { addTaskToBrowserMainThread, copyTextToClipboard } from '@swanlab-vue/utils/browser'
import { formatNumber2SN } from '@swanlab-vue/utils/common'
import { t } from '@swanlab-vue/i18n'
import { getTimes } from '@swanlab-vue/utils/time'
import { isApple } from '@swanlab-vue/utils/browser'
import { message } from '@swanlab-vue/components/message'
import makeSmooth from './smooth'
import { needSmooth } from './smooth'
import LineChartTooltip from '../components/LineChartTooltip.vue'
import LineChartLegend from '../components/LineChartLegend.vue'
import '@swanlab-vue/docs'

// ---------------------------------- 配置 ----------------------------------
const props = defineProps({
  title: {
    type: String,
    required: true
  },
  chart: {
    type: Object,
    required: true
  },
  index: {
    type: Number,
    required: true
  }
})

// ---------------------------------- 错误处理，如果chart.error存在，则下面的api都将不应该被执行 ----------------------------------
// 数据源 array
const source = props.chart.source
// source的长度如果等于error的长度，说明所有数据都有问题,取第一个的error即可
const error = ref(
  props.chart.error && source.length === Object.keys(props.chart.error).length ? props.chart.error[source[0]] : null
)
// 图表模式，multi或者single
const multi = props.chart.multi
// 当前命名空间下除了自己的所有chartRef
const chartRefList = inject('chartRefList')
const chartRefListExceptSelf = computed(() => {
  // 将列表中除了props.index的所有chartRef过滤出来
  return chartRefList.value.filter((_, i) => i !== props.index)
})
// 平滑配置
let smoothMethod = inject('smoothMethod')?.value || undefined
// ---------------------------------- 图表样式配置 ----------------------------------
// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const getColor = inject('getColor')
const defaultColor = inject('defaultColor')
const rootStyle = getComputedStyle(document.documentElement)
// 边框颜色，通过js获取css变量值
const borderColor = rootStyle.getPropertyValue('--outline-default')
// 网格线颜色，通过js获取css变量值
const gridColor = rootStyle.getPropertyValue('--outline-dimmest')
// 十字准线颜色，通过js获取css变量值
const crosshairsColor = rootStyle.getPropertyValue('--primary-dimmest')
// 线段默认宽度
const lineWidth = 1.5
// 线段加粗宽度
const thickerLineWidth = 3.5

// ---------------------------------- 样式注册，数据点样式注册，如果是最后一个，会放大 ----------------------------------
G2.registerShape('point', 'last-point', {
  draw(cfg, container) {
    const point = { x: cfg.x, y: cfg.y }
    // console.log('point', cfg.data)
    const shape = container.addShape('circle', {
      name: 'point',
      attrs: {
        x: point.x,
        y: point.y,
        fill: cfg.color || 'red',
        opacity: cfg?.data?._last ? 1 : 0,
        r: 3
      }
    })
    return shape
  }
})

// ---------------------------------- 组件渲染逻辑 ----------------------------------
// 组件对象
const lineChartRef = ref()
const g2Ref = ref()
const g2ZoomRef = ref()
const tooltipRef = ref()
const tooltipZoomRef = ref()
// 参考字段和显示名称
const reference = props.chart.reference
// 拿到参考系，未来图表可能有不同的x轴依据，比如step、time等，这里需要根据设置的reference来决定
const { xField, xTitle } = UTILS.refrence2XField[reference]
// 默认y轴的依据key是data
const yField = 'data'
const seriesField = 'series'
const colorField = 'type'
/**
 * 创建图表函数
 * @param { HTMLElement } dom 图表挂载的dom
 * @param { Array } data 图表数据
 * @param { Object } config 图表的一些其他配置
 * @param { bool } zoom 是否放大
 */
const createChart = (dom, data, config = {}, zoom = false) => {
  if (!dom) return
  const c = new Line(dom, {
    data,
    // 默认的x轴依据key为step
    xField,
    // 默认的y轴依据key为data
    yField,
    // 多数据的时候，需要设置seriesField，单数据也可以设置，但是不希望出现label
    // seriesField,
    colorField,
    // 自己写图例
    legend: false,
    // 多数据的时候颜色通过回调拿到
    color: ({ series }) => {
      // 是否有&smooth后缀
      const smooth = series.includes('&smooth')
      series = series.replace('&smooth', '')
      // 如果没有&smooth后缀，直接返回颜色
      const color = getColor(series, source.indexOf(series), props.chart.id)
      if (!smooth) return color
      else return UTILS.transparentColor(color)
    },
    point: {
      shape: 'last-point'
    },
    lineStyle: {
      lineWidth
    },
    // 坐标轴相关
    xAxis: {
      // 自定义坐标轴的刻度，暂时没有找到文档，通过源码来看是返回一个数组，数组内是字符串，代表刻度
      tickCount: 5,
      type: 'linear',
      // 在此处完成X轴数据的格式化
      label: {
        formatter: (data) => {
          // console.log('data', data)
          // 如果是小数，返回空
          if (data % 1 !== 0) return ''
          // 如果是100的倍数且大于1000，返回k
          if (data % 100 === 0 && data >= 1000) return `${data / 1000}k`
          return data
        }
      },

      // x轴坐标轴样式
      line: {
        style: {
          stroke: borderColor,
          lineWidth: 2
        }
      },
      // x轴刻度样式
      tickLine: {
        length: 4,
        style: {
          stroke: borderColor,
          lineWidth: 2
        }
      }
    },
    yAxis: {
      min: null,
      label: {
        // 在此处完成Y轴数据的格式化
        formatter: (data) => {
          return formatNumber2SN(data)
        }
      },
      // y轴坐标轴样式
      line: {
        style: {
          stroke: borderColor,
          lineWidth: 2
        }
      },
      // y轴刻度样式
      tickLine: {
        length: 4,
        style: {
          stroke: borderColor,
          lineWidth: 2
        }
      },
      // 网格线
      grid: {
        line: {
          style: {
            stroke: gridColor
          }
        }
      }
    },
    // 悬浮提示相关
    tooltip: {
      // 在此处完成悬浮数据提示的格式化
      // 如果需要自定义浮窗，可以用下面的customContent
      // formatter: (data) => {
      //   // console.log(data)
      //   // 如果data.series是undefined，说明是单数据,直接显示source[0]即可
      //   const name = data.series ? data.series : source[0]
      //   return { name, value: formatNumber2SN(data.data) }
      // },
      follow: true,
      enterable: false,
      shared: true,
      position: 'top',
      showMarkers: true,
      customContent: () => '',
      domStyles: {
        // 自己写tooltip的样式
        'g2-tooltip': {
          boxShadow: 'none',
          borderWidth: 'none',
          borderRadius: 'none'
        }
      },
      showCrosshairs: true,
      crosshairs: {
        line: {
          style: {
            stroke: crosshairsColor,
            lineWidth: 2
          }
        }
      }
    },
    // 大小相关
    height: originalHeight,
    width: undefined,
    autoFit: true,
    // 开启一些交互
    // interactions: [{ type: 'element-active' }],
    // 平滑曲线
    smooth: false,
    animation: false,
    ...config
  })
  c.render()
  addTaskToBrowserMainThread(() => {
    registerTooltipEvent(dom, zoom)
  })
  // 监听鼠标移入事件
  // c.on('plot:mouseenter', (evt) => {
  // })
  // 监听鼠标移出事件
  c.on('plot:mouseleave', () => {
    restoreByTag(c, zoom, nowThickenTag, nowThickenColor)
  })
  // 监听鼠标移动事件
  c.on('plot:mousemove', (evt) => {
    // console.log('plot:mousemove', evt)
    if (!nowData) return
    const y = evt.y
    // 遍历所有的nowData，寻找其中与y绝对值最小的点
    let index = undefined
    let min = Infinity
    for (let i = 0; i < nowData.length; i++) {
      const data = nowData[i]
      const d = Math.abs(data.y - y)
      if (d < min && !nowData[i].data.series.includes('&smooth')) {
        min = d
        index = i
      }
    }
    // console.log('鼠标移动', props.chart.name)
    if (index !== undefined) thickenByTag(c, zoom, nowData[index].data.series, nowData[index].color)
  })

  return c
}

// ---------------------------------- 数据格式化 ----------------------------------
/**
 * 为了将数据格式化为图表可用的格式，需要将数据源中的数据进行格式化
 * 遍历data的所有key，合并其中的list为一个数组
 * @param { OpenChartData } data 待格式化的数据
 * @returns { Object } 格式化后的数据, { d: [{}, {}, ...], config: {} } config是图表的一些其他配置
 */
const format = (data) => {
  // 如果配置了平滑方法，需要动态计算平滑后的数据
  if (smoothMethod) {
    const smoothData = {}
    for (const key in data) {
      if (!data[key]) continue
      // data中出现了字符串
      console.log('data[key]', data[key])
      if (data[key].list.some((item) => typeof item.data === 'string')) {
        smoothData[key] = data[key]
        continue
      }
      smoothData[key] = {
        ...data[key],
        list: makeSmooth(data[key].list, smoothMethod)
      }
    }
    // 将原本的data的key添加&smooth后，但是注意不要修改原本的data，应该新建一个对象
    for (const key in data) {
      smoothData[`${key}&smooth`] = data[key]
    }
    data = smoothData
  }
  // console.log('data', data)
  // 如果source的长度小于1，抛出错误
  // console.log(data)
  if (source.length < 1) throw new Error('source length must be greater than 1')
  // 新的数据,遍历得到
  const d = []
  Object.keys(data).forEach((key) => {
    // 如果key存在于props.chart.error中，说明这个数据有问题，直接返回
    if (props.chart.error && props.chart.error[key]) return
    // 如果不是单数据，需要将所有数据的list合并为一个数组
    if (!data[key]) return
    data[key].list.forEach((item) => {
      // item新加series字段，用于标识数据来源
      d.push({ ...item, series: key })
    })
  })
  // console.log('d', d)
  // 依据xField排序，从小到大
  d.sort((a, b) => a[xField] - b[xField])
  // console.log('d', d)
  // console.log('data', data)
  // 依据source生成legend，排序依据是source，数据结构：[{name, color}]
  const items = []
  for (const s of source) {
    if (props.chart.error && props.chart.error[s]) continue
    if (!data[s]) continue
    items.push({ name: s, color: getColor(s, source.indexOf(s), props.chart.id), experiment_id: data[s].experiment_id })
  }
  legend.value = items
  // console.log('legend', legend.value)
  return { d, config: multi || smoothMethod ? { seriesField } : { color: defaultColor } }
}

const formatTime = (time) => {
  let { year, month, day, hour, minute } = getTimes(time)
  // year只保留后两位
  year = year.toString().slice(-2)
  // 格式化minute，如果是个位数，前面加0
  minute = minute.toString().length === 1 ? `0${minute}` : minute
  let pm = false
  if (hour > 12) {
    hour -= 12
    pm = true
  }
  let result = `${year}/${month}/${day} ${hour}:${minute}`
  if (pm) result += ' PM'
  else result += ' AM'
  return result
}

provide('formatTime', formatTime)
provide('formatNumber2SN', formatNumber2SN)
// ---------------------------------- 渲染、重渲染功能 ----------------------------------
let chartObj = null
// 渲染
const render = (data) => {
  // console.log('渲染折线图')
  // console.log('data', data)
  const { d, config } = format(data)
  // console.log('data', data)
  chartObj && chartObj.destroy()
  chartObj = createChart(g2Ref.value, d, config)
  // console.log('chartObj', chartObj)
  // 可以使用update api来更新配置
}
// 重渲染
const change = (data) => {
  const { d } = format(data)
  chartObj.changeData(d)
  // // change函数等于render函数
  // chartObj.destroy()
  // chartObj = createChart(g2Ref.value, d, { animation: false, ...config })
}

// ---------------------------------- 放大功能 ----------------------------------
// 是否放大
const isZoom = ref(false)
let zoomChartObj = null
// 放大数据
const zoom = (data) => {
  isZoom.value = true
  // 当前window的高度
  const { d, config } = format(data)
  const height = window.innerHeight * 0.6
  addTaskToBrowserMainThread(() => {
    zoomChartObj = createChart(
      g2ZoomRef.value,
      d,
      {
        interactions: [{ type: 'brush-x' }],
        height,
        ...config
      },
      true
    )
  })
}

// ---------------------------------- 多linechart之间的tooltip联动 ----------------------------------

let manual = true
// 调用此方法则必然是自动触发
const lineShowTooltip = (point) => {
  if (error.value) return
  manual = false
  // console.log('lineShowTooltip', props.index)
  // console.log('title', title)
  chartObj?.chart?.showTooltip(point)
}
const lineHideTooltip = () => {
  if (error.value) return
  manual = false
  chartObj?.chart?.hideTooltip()
}
let point = null
const registerTooltipEvent = (dom, zoom) => {
  // 给 tooltip 添加点击事件
  if (error.value) return
  // 代理实现tooltip的显示和隐藏
  const chart = zoom ? zoomChartObj : chartObj
  const tRef = zoom ? tooltipZoomRef : tooltipRef

  chart.on('tooltip:show', (evt) => {
    // console.log('tooltip:show', evt)
    // 获取当前tooltip的left
    const items = evt.data.items
    addTaskToBrowserMainThread(() => {
      const left = dom.querySelector('.g2-tooltip').style.left
      tRef.value.show(items, dom.clientWidth, left)
    })
    // 非js触发的tooltip，不执行下面的逻辑
    if (!manual) {
      // console.log('auto show tooltip')
      return
    }
    nowData = evt.data.items
    // console.log('nowData', nowData)
    point = { x: evt.data.x, y: evt.data.y }
    // 通知其他图表，当前图表的数据被hover到了
    !zoom &&
      chartRefListExceptSelf.value?.forEach((chart) => {
        chart?.showTooltip(point)
      })
    manual = true
  })
  chart.on('tooltip:hide', (...args) => {
    // console.log('tooltip:hide', chart)
    // 非js触发的tooltip，不执行下面的逻辑
    tRef.value.hide()
    if (manual && !zoom)
      // 通知其他图表，当前图表的数据被hover到了
      chartRefListExceptSelf.value?.forEach((chart) => {
        chart?.hideTooltip(...args)
      })
    nowData = null
    manual = true
    point = null
  })
  // brush相关事件
  chart.on(G2.BRUSH_FILTER_EVENTS.BEFORE_FILTER, (...agrs) => {
    console.log('BEFORE_FILTER', agrs)
  })
}

// ---------------------------------- tooltip出现并且是非js触发时，可执行copy操作 ----------------------------------
// 当前tooltip的数据,用于copy
let nowData = null
// 全局注册keydown事件，当mac端触发command+c，windows端触发ctrl+c时，且nowData不为null，执行copy操作
const handleCopy = (e) => {
  if (error.value) return
  if (nowData === null) {
    return
    // return console.log('非js触发的tooltip，或者未悬浮，不执行copy操作')
  }
  if (e.key === 'c' && (isApple ? e.metaKey : e.ctrlKey)) {
    e.preventDefault()
    // console.log('copy:', nowData)
    // nowData依据data降序
    nowData.sort((a, b) => b.data.data - a.data.data)
    // 生成copy的内容，zoom和非zoom样式不一样
    let content = ''
    if (!isZoom.value) {
      // 复制数据和tag：tag data
      for (const d of nowData) {
        content += `${d.data.series} ${formatNumber2SN(d.data.data)}\n`
      }
    } else {
      // 复制数据、时间和tag：tag data time
      for (const d of nowData) {
        content += `${d.data.series} ${formatNumber2SN(d.data.data)} ${formatTime(d.data.create_time)}\n`
      }
    }
    copyTextToClipboard(content, () => message.success(t('chart.charts.line.copy.success')))
  }
}
window.addEventListener('keydown', handleCopy)
onUnmounted(() => {
  window.removeEventListener('keydown', handleCopy)
})

// ---------------------------------- 控制线段加粗 ----------------------------------
/**
 * 当前加粗的tag,需要注意的是当前加粗的tag可能不存在于当前图表的数据中
 * multi为true时，按照tag加粗，否则按照color加粗
 */
let nowThickenTag = null
/**
 * 当前加粗的颜色，需要注意的是当前加粗的颜色可能不存在于当前图表的数据中
 * mutlti为false时，按照color加粗，否则按照tag加粗
 */
let nowThickenColor = null
/**
 * 加粗指定tag的线段
 * @param { Object } plot plot对象，chartObj或者zoomChartObj
 * @param { bool } zoom 是否放大
 * @param { string } tag 实验标签
 */
const thickenByTag = (plot, zoom, tag, color) => {
  // console.log('触发:', props.chart.name)
  if (!tag || !color) return
  if (!plot) return
  if (tag === nowThickenTag) return // 避免重复加粗，会造成无限回调
  if (tag !== nowThickenTag && nowThickenTag) restoreByTag(plot, zoom, nowThickenTag, nowThickenColor)
  const els = plot.chart.getElements()
  for (const e of els) {
    // 多数据模式下，依据tag加粗
    if (multi) {
      // 如果是array，取[0],如果是object，取series
      const series = e.model.data[0]?.series || e.model.data.series
      // console.log('series', series)
      if (series === tag) {
        // console.log('加粗', props.chart.name, e)
        e.update({ ...e.model, style: { lineWidth: thickerLineWidth } })
        break
      }
    }
    // 单数据模式下，依据color加粗
    else {
      if (color === e.model.color) {
        e.update({ ...e.model, style: { lineWidth: thickerLineWidth } })
        break
      }
    }
  }
  nowThickenTag = tag
  nowThickenColor = color
  // 加粗其他图表的数据，保持联动
  chartRefListExceptSelf.value.forEach((chart) => {
    chart?.thicken(zoom, tag, color)
  })
}

/**
 * 将指定颜色的线段恢复原状
 * @param { Object } plot plot对象，chartObj或者zoomChartObj
 * @param { bool } zoom 是否放大
 * @param { string } tag 实验标签
 */
const restoreByTag = (plot, zoom, tag, color) => {
  if (!tag || !color) return
  if (!nowThickenTag) return
  const els = plot.chart.getElements()
  for (const e of els) {
    // 多数据模式下，依据tag恢复
    if (multi) {
      const series = e.model.data[0]?.series || e.model.data.series
      if (series === tag) {
        e.update({ ...e.model, style: { lineWidth: lineWidth } })
        break
      }
    }
    // 单数据模式下，依据color恢复
    else {
      if (color === e.model.color) {
        e.update({ ...e.model, style: { lineWidth: lineWidth } })
        break
      }
    }
  }
  nowThickenTag = null
  nowThickenColor = null
  chartRefListExceptSelf.value.forEach((chart) => {
    chart?.thin(zoom, tag, color)
  })
}

/**
 * 用于交给其他图表联动调用来加粗线段
 * @param { Object } plot plot对象，chartObj或者zoomChartObj
 * @param { bool } zoom 是否放大
 * @param { string } color 颜色
 */
const thickenByTagLinkage = (zoom, tag, color) => {
  const plot = zoom ? zoomChartObj : chartObj
  thickenByTag(plot, zoom, tag, color)
}

/**
 * 用于交给其他图表联动调用来恢复线段粗细
 */
const restoreByTagLinkage = (zoom, tag, color) => {
  const plot = zoom ? zoomChartObj : chartObj
  restoreByTag(plot, zoom, tag, color)
}

// ---------------------------------- 图例渲染，在这保存数据，传给legend组件 ----------------------------------

const legend = ref(null)

const legendHoverin = (item, zoom) => {
  // console.log('mouseenter', item)
  const plot = zoom ? zoomChartObj : chartObj
  thickenByTag(plot, zoom, item.name, item.color)
}

const legendHoverout = (item, zoom) => {
  // console.log('mouseleave', item)
  const plot = zoom ? zoomChartObj : chartObj
  restoreByTag(plot, zoom, item.name, item.color)
}

// ---------------------------------- 图表平滑 ----------------------------------
/**
 * @type {OpenChartData} chartData
 */
const chartData = inject('data')
const smooth = (method) => {
  console.log('smooth by linechart:', chartData)
  // console.log('smooth by linechart:', method)
  if (!needSmooth(method)) {
    smoothMethod = null
  } else {
    smoothMethod = method
  }
  render(chartData)
}

// ---------------------------------- 动态设置原始高度 ----------------------------------

let originalHeight = 220
/**
 * @type {setOriginalChartHeight} setLineChartHeight
 */
const setLineChartHeight = (height, maxHeight = 800, minHeight = 200) => {
  if (height > maxHeight) return
  if (height < minHeight) return
  const nowLineChartHeight = lineChartRef.value.clientHeight
  const nowLineChartCanvas = lineChartRef.value.getElementsByTagName('canvas')[0]
  if (!nowLineChartCanvas) return (originalHeight = height - nowLineChartCanvasHeight)
  const nowLineChartCanvasHeight = nowLineChartCanvas.clientHeight
  console.log('nowLineChartHeight', nowLineChartHeight)
  console.log('nowLineChartCanvasHeight', nowLineChartCanvasHeight)
  const staticHeight = nowLineChartHeight - nowLineChartCanvasHeight
  // 设置canvas高度，需要减去静态高度
  nowLineChartCanvas.style.height = `${height - staticHeight}px`
}

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  render,
  change,
  zoom,
  smooth,
  showTooltip: lineShowTooltip,
  hideTooltip: lineHideTooltip,
  thicken: thickenByTagLinkage,
  thin: restoreByTagLinkage,
  setHeight: setLineChartHeight
})
</script>

<style lang="scss"></style>
