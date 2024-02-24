<template>
  <!-- 图表标题 -->
  <p class="text-center font-semibold select-none">{{ title }}</p>
  <div class="flex flex-col justify-center grow text-dimmer gap-2" v-if="error">
    <SLIcon class="mx-auto h-5 w-5" icon="error" />
    <p class="text-center text-sm">
      {{ $t('common.chart.charts.line.error', { type: error['data_class'], tag: source[0] }) }}
    </p>
  </div>
  <template v-else>
    <!-- x轴坐标单位 -->
    <p class="absolute right-5 bottom-10 text-xs text-dimmer scale-90 select-none">{{ xTitle }}</p>
    <!-- 图表主体 -->
    <div ref="g2Ref"></div>
    <!-- 放大效果 -->
    <SLModal class="p-10 pt-0 overflow-hidden" max-w="-1" v-model="isZoom">
      <p class="text-center mt-4 mb-10 text-2xl font-semibold">{{ title }}</p>
      <div ref="g2ZoomRef"></div>
      <p class="absolute right-12 bottom-16 text-xs text-dimmer scale-90 select-none">{{ xTitle }}</p>
    </SLModal>
  </template>
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
import { ref, inject, computed, onUnmounted } from 'vue'
import { addTaskToBrowserMainThread, copyTextToClipboard } from '@swanlab-vue/utils/browser'
import { formatNumber2SN } from '@swanlab-vue/utils/common'
import { t } from '@swanlab-vue/i18n'
import { getTimes } from '@swanlab-vue/utils/time'
import { isApple } from '@swanlab-vue/utils/browser'
import { message } from '@swanlab-vue/components/message'

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
// 数据源 arrya
const source = props.chart.source
// source的长度如果等于error的长度，说明所有数据都有问题,取第一个的error即可
const error = ref(source.length === Object.keys(props.chart.error).length ? props.chart.error[source[0]] : null)
// 图表模式，mutli或者single
const mutli = props.chart.mutli
// ---------------------------------- 图表样式配置 ----------------------------------
// 后续需要适配不同的颜色，但是Line不支持css变量，考虑自定义主题或者js获取css变量完成计算
const colors = inject('colors')
if (!colors) throw new Error('colors is not defined, please provide colors in parent component')
if (!colors.getSeriesColor) throw new Error('colors.getSeries is not defined, please provide getSeries in colors')
const rootStyle = getComputedStyle(document.documentElement)
// 边框颜色，通过js获取css变量值
const borderColor = rootStyle.getPropertyValue('--outline-default')
// 网格线颜色，通过js获取css变量值
const gridColor = rootStyle.getPropertyValue('--outline-dimmest')
// tooltip内容，从左到右，分别是颜色，步长，值，时间，tag，先第一行
let tooltipContent = `<div class="lc-tooltip-item-zoom">`
for (const c of [
  ['', 'lc-tooltip-color'],
  [t('common.chart.charts.share.step'), 'lc-tooltip-step'],
  [t('common.chart.charts.share.value'), 'lc-tooltip-value'],
  [t('common.chart.charts.share.time'), 'lc-tooltip-time'],
  [t('common.chart.charts.share.tag'), 'lc-tooltip-tag']
]) {
  tooltipContent += `<p class="${c[1]}">${c[0]}</p>`
}
tooltipContent += `</div>`
// tooltip的x轴偏移量
const tooltipXOffset = 50

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
const g2Ref = ref()
const g2ZoomRef = ref()
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
  const c = new Line(dom, {
    data,
    // 默认的x轴依据key为step
    xField,
    // 默认的y轴依据key为data
    yField,
    // 多数据的时候，需要设置seriesField，单数据也可以设置，但是不希望出现label
    // seriesField,
    colorField,
    legend: {
      // flipPage: false,
      pageNavigator: {
        marker: {
          style: {
            // 非激活，不可点击态时的填充色设置
            inactiveFill: '#000',
            inactiveOpacity: 0.45,
            // 默认填充色设置
            fill: '#000',
            opacity: 0.8,
            size: 8
          }
        },
        text: {
          style: {
            fill: '#ccc',
            fontSize: 8
          }
        }
      }
    },
    // 多数据的时候颜色通过回调拿到，colors应该自带getSeries方法
    color: ({ series }) => {
      return colors.getSeriesColor(series, source.indexOf(series))
    },
    point: {
      shape: 'last-point'
    },
    // 坐标轴相关
    xAxis: {
      // 自定义坐标轴的刻度，暂时没有找到文档，通过源码来看是返回一个数组，数组内是字符串，代表刻度
      tickCount: 5,
      type: 'linear',
      // tickMethod: (cfg) => {
      //   // console.log('chart', props.chart.name)
      //   console.log('cfg', cfg)
      //   console.log('cfg.values', cfg.values)
      //   // 数据个数
      //   const num = cfg.values?.length || 0
      //   // 每隔多少个显示一个刻度
      //   const step = Math.ceil(num / 5)
      //   console.log('step', step)
      //   const tick = []
      //   let i = 0
      //   while (i < num) {
      //     tick.push(i)
      //     i += step
      //   }
      //   console.log('tick', tick)
      //   return tick
      // },
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
    // 图例相关
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
      position: 'left',
      offset: 10,
      customContent: (title, data) => {
        // 网格布局，一行五列
        // console.log('data', data)
        // 一行一行生成,放大样式和非放大样式不一样
        let content = zoom ? tooltipContent : ``
        // data 根据value降序
        data.sort((a, b) => b.data.data - a.data.data)
        addTaskToBrowserMainThread(() => {
          const g2Tooltip = dom.querySelector('.g2-tooltip')
          console.log('g2Tooltip top', g2Tooltip.style.top)
          const tooltipDom = g2Ref.value.querySelector('.lc-tooltip')
          tooltipDom.style.top = '-' + g2Tooltip.style.top
          // 设置left或者right的值，这将取决于tooltip左右两侧是否超出dom的宽度
          // 如果左侧超出，设置right，否则设置left
          const clientWidth = parseFloat(g2Ref.value.clientWidth)
          const left = parseFloat(g2Tooltip.style.left)
          console.log('clientWidth', clientWidth)
          console.log('left', left)
          if (left + tooltipXOffset + tooltipDom.clientWidth > clientWidth) {
            tooltipDom.style.right = `${tooltipXOffset - 50}px`
          } else {
            tooltipDom.style.left = `${tooltipXOffset}px`
          }
        })
        if (!zoom) {
          // 缩小样式
          for (const d of data) {
            content += `<div class="lc-tooltip-item-no-zoom" style="color:${d.color}">`
            // 先颜色
            content += `<span class="lc-tooltip-color lc-tooltip-color-rect"></span>`
            // 再步长
            content += `<span class="lc-tooltip-step">${d.data.index}</span>`
            // 再值
            content += `<span class="lc-tooltip-value">${formatNumber2SN(d.data.data)}</span>`
            // 再tag
            content += `<span class="lc-tooltip-tag">${d.data.series}</span>`
            content += `</div>`
          }
        } else {
          for (const d of data) {
            content += `<div class="lc-tooltip-item-zoom" style="color:${d.color}">`
            // 先颜色
            content += `<span class="lc-tooltip-color lc-tooltip-color-rect"></span>`
            // 再步长
            content += `<span class="lc-tooltip-step">${d.data.index}</span>`
            // 再值
            content += `<span class="lc-tooltip-value">${formatNumber2SN(d.data.data)}</span>`
            // 再时间
            content += `<span class="lc-tooltip-time">${formatTime(d.data.create_time)}</span>`
            // 再tag
            content += `<span class="lc-tooltip-tag">${d.data.series}</span>`
            content += `</div>`
          }
        }
        // 最后加一行tip
        content += `<p class="lc-tooltip-tip">${tip}</p>`
        return `<div class="lc-tooltip">${content}</div>`
      },
      domStyles: {
        'g2-tooltip': {
          boxShadow: 'none',
          borderWidth: 'none',
          borderRadius: 'none'
        }
      }
    },
    // 大小相关
    height: 220,
    width: undefined,
    autoFit: true,
    // 开启一些交互
    interactions: [{ type: 'hover-cursor' }],
    // 平滑曲线
    smooth: false,
    ...config
  })
  c.render()
  addTaskToBrowserMainThread(() => {
    registerTooltipEvent(zoom)
  })
  // 监听鼠标移入事件
  c.on('plot:mouseenter', (evt) => {
    // console.log('plot:mouseenter', evt)
    // console.log(c.chart.getElements())
    // const e = c.chart.getElements()[0]
    // c.chart.getElements()[0].update({ ...e.model, defaultStyle: { lineWidth: 5 } })
    // 通过颜色判断类型，只需要改变每个类型的第一个元素, color => e
    const eMap = new Map()
    c.chart.getElements().forEach((e) => {
      if (!eMap.has(e.model.color)) eMap.set(e.model.color, e)
    })
    eMap.forEach((e) => {
      e.update({ ...e.model, defaultStyle: { lineWidth: 5 } })
    })
  })
  // 监听鼠标移出事件
  c.on('plot:mouseleave', (evt) => {
    // console.log('plot:mouseleave', evt)
    // 同理
    // const eMap = new Map()
    // c.chart.getElements().forEach((e) => {
    //   if (!eMap.has(e.model.color)) eMap.set(e.model.color, e)
    // })
    // eMap.forEach((e) => {
    //   e.update({ ...e.model, defaultStyle: { lineWidth: 2 } })
    // })
  })

  return c
}

// ---------------------------------- 数据格式化 ----------------------------------
/**
 * 为了将数据格式化为图表可用的格式，需要将数据源中的数据进行格式化
 * 遍历data的所有key，合并其中的list为一个数组
 * @param { Object } data 待格式化的数据
 * @returns { Object } 格式化后的数据, { d: [{}, {}, ...], config: {} } config是图表的一些其他配置
 */
const format = (data) => {
  // 如果source的长度小于1，抛出错误
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
  // 如果source的长度大于1，需要设置seriesField
  return { d, config: mutli ? { seriesField } : { color: colors[0] } }
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
// ---------------------------------- 渲染、重渲染功能 ----------------------------------
let chartObj = null
// 渲染
const render = (data) => {
  // console.log('渲染折线图')
  // console.log('data', data)
  const { d, config } = format(data)
  // console.log('data', data)
  chartObj = createChart(g2Ref.value, d, config)
  // console.log('chartObj', chartObj)
  // 可以使用update api来更新配置
}
// 重渲染
const change = (data) => {
  const { d, config } = format(data)
  // change函数等于render函数
  chartObj.destroy()
  chartObj = createChart(g2Ref.value, d, { animation: false, ...config })
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

// ---------------------------------- 额外功能：多linechart之间的tooltip联动 ----------------------------------
const chartsRefList = inject('chartsRefList')
const lineChartsRef = computed(() => {
  // 将列表中除了props.index的所有chartRef过滤出来
  return chartsRefList.value.filter((item, i) => i !== props.index && item.chartRef?.lineShowTooltip)
})

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

const registerTooltipEvent = (zoom) => {
  // 给 tooltip 添加点击事件
  if (error.value) return
  const chart = zoom ? zoomChartObj : chartObj
  chart.on('tooltip:show', (evt) => {
    // 非js触发的tooltip，不执行下面的逻辑
    if (!manual) {
      // console.log('auto show tooltip')
      return
    }
    nowData = evt.data.items
    // console.log('nowData', nowData)
    const point = { x: evt.data.x, y: evt.data.y }
    // 通知其他图表，当前图表的数据被hover到了
    !zoom &&
      lineChartsRef.value?.forEach((chart) => {
        chart.chartRef.lineShowTooltip(point)
      })
    manual = true
  })
  chart.on('tooltip:hide', (...args) => {
    // 非js触发的tooltip，不执行下面的逻辑
    if (manual && !zoom)
      // 通知其他图表，当前图表的数据被hover到了
      lineChartsRef.value?.forEach((chart) => {
        chart.chartRef.lineHideTooltip(...args)
      })
    nowData = null
    manual = true
  })
}

// ---------------------------------- tooltip出现并且是非js触发时，可执行copy操作 ----------------------------------
const tip = isApple ? t('common.chart.charts.line.copy.apple') : t('common.chart.charts.line.copy.windows')
// 当前tooltip的数据,用于copy
let nowData = null
// 全局注册keydown事件，当mac端触发command+c，windows端触发ctrl+c时，且nowData不为null，执行copy操作
const handelCopy = (e) => {
  if (error.value) return
  if (nowData === null) {
    // return
    return console.log('非js触发的tooltip，或者未悬浮，不执行copy操作')
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
    copyTextToClipboard(content, () => message.success(t('common.chart.charts.line.copy.success')))
  }
}
window.addEventListener('keydown', handelCopy)
onUnmounted(() => {
  window.removeEventListener('keydown', handelCopy)
})

// ---------------------------------- 暴露api ----------------------------------
defineExpose({
  render,
  change,
  zoom,
  lineShowTooltip,
  lineHideTooltip
})
</script>

<style lang="scss">
.lc-tooltip {
  @apply py-2 px-3 absolute bg-default border rounded;
  box-shadow: rgba(21, 24, 31, 0.16) 0px 12px 24px 0px;
  p {
    @apply text-xs text-default font-semibold;
  }
  .lc-tooltip-item-no-zoom,
  .lc-tooltip-item-zoom {
    @apply flex items-center gap-3;
    &:not(:last-child) {
      @apply mb-1.5;
    }
    .lc-tooltip-color {
      @apply w-5 flex items-center;
    }
    .lc-tooltip-color-rect {
      &::before {
        content: '';
        display: inline-block;
        width: 20px;
        height: 6px;
        border-radius: 2px;
        margin-right: 5px;
        background-color: currentColor;
      }
    }
  }
  .lc-tooltip-tip {
    @apply font-normal text-dimmest text-xs;
  }
}

.lc-tooltip-item-no-zoom {
  .lc-tooltip-step {
    @apply font-semibold;
    &::after {
      content: ':';
      @apply font-semibold;
    }
  }
  .lc-tooltip-value {
    @apply w-10 text-left font-semibold;
  }
  .lc-tooltip-tag {
    @apply truncate;
    max-width: 128px;
  }
}

.lc-tooltip-item-zoom {
  .lc-tooltip-step {
    @apply w-7;
  }
  .lc-tooltip-value {
    @apply col-span-1;
    @apply w-10 text-left;
  }
  .lc-tooltip-time {
    @apply w-28;
  }

  .lc-tooltip-tag {
    @apply truncate;
    max-width: 160px;
  }
}
</style>
