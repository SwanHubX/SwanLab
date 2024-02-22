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
import { ref, inject, computed } from 'vue'
import { addTaskToBrowserMainThread } from '@swanlab-vue/utils/browser'
import { formatNumber2SN } from '@swanlab-vue/utils/common'
import { t } from '@swanlab-vue/i18n'
import { getTimes } from '@swanlab-vue/utils/time'

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
let tooltipContent = ``
for (const c of [
  ['', 'lc-tooltip-color'],
  [t('common.chart.charts.share.step'), 'lc-tooltip-step'],
  [t('common.chart.charts.share.value'), 'lc-tooltip-value'],
  [t('common.chart.charts.share.time'), 'lc-tooltip-time'],
  [t('common.chart.charts.share.tag'), 'lc-tooltip-tag']
]) {
  tooltipContent += `<p class="${c[1]}">${c[0]}</p>`
}

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
// 创建图表的函数
const createChart = (dom, data, config = {}) => {
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
      // 如果需要自定义浮窗，可以用下面的customContent，但是目前不管
      // formatter: (data) => {
      //   // console.log(data)
      //   // 如果data.series是undefined，说明是单数据,直接显示source[0]即可
      //   const name = data.series ? data.series : source[0]
      //   return { name, value: formatNumber2SN(data.data) }
      // }
      customContent: (title, data) => {
        // 网格布局，一行五列
        console.log('data', data)
        // 一行一行生成
        let content = tooltipContent
        for (const d of data) {
          // 先颜色
          content += `<span class="lc-tooltip-color lc-tooltip-color-rect" style="color:${d.color}"></span>`
          // 再步长
          content += `<span class="lc-tooltip-step" style="color:${d.color}">${d.data.index}</span>`
          // 再值
          content += `<span class="lc-tooltip-value" style="color:${d.color}">${formatNumber2SN(d.data.data)}</span>`
          // 再时间
          content += `<span class="lc-tooltip-time" style="color:${d.color}">${formatTime(d.data.create_time)}</span>`
          // 再tag
          content += `<span class="lc-tooltip-tag" style="color:${d.color}">${d.data.series}</span>`
        }
        return `<div class="lc-tooltip">${content}</div>`
      }
    },
    // 大小相关
    height: 220,
    width: undefined,
    autoFit: true,
    // 开启一些交互
    interactions: [{ type: 'element-active' }],
    // 平滑曲线
    smooth: false,
    ...config
  })
  c.render()
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
  console.log('d', d)
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
  let result = `${month}/${day}/${year} ${hour}:${minute}`
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
  registerTooltipEvent()
}
// 重渲染
const change = (data) => {
  const { d, config } = format(data)
  // change函数等于render函数
  chartObj.destroy()
  chartObj = createChart(g2Ref.value, d, { animation: false, ...config })
  registerTooltipEvent()
}

// ---------------------------------- 放大功能 ----------------------------------
// 是否放大
const isZoom = ref(false)
// 放大数据
const zoom = (data) => {
  isZoom.value = true
  // 当前window的高度
  const { d, config } = format(data)
  const height = window.innerHeight * 0.6
  addTaskToBrowserMainThread(() => {
    createChart(g2ZoomRef.value, d, {
      interactions: [{ type: 'brush-x' }, { type: 'element-active' }],
      height,
      ...config
    })
  })
}

// ---------------------------------- 额外功能 ----------------------------------
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

const registerTooltipEvent = () => {
  // 给 tooltip 添加点击事件
  if (error.value) return
  chartObj.on('tooltip:show', (evt) => {
    if (!manual) {
      // console.log('auto show tooltip')
      return
    }
    const point = { x: evt.data.x, y: evt.data.y }
    // 通知其他图表，当前图表的数据被hover到了
    lineChartsRef.value?.forEach((chart) => {
      chart.chartRef.lineShowTooltip(point)
    })
    manual = true
  })
  chartObj.on('tooltip:hide', (...args) => {
    // 通知其他图表，当前图表的数据被hover到了
    lineChartsRef.value?.forEach((chart) => {
      if (!manual) {
        return
        // return console.log('auto hide tooltip')
      }
      chart.chartRef.lineHideTooltip(...args)
    })
    manual = true
  })
}
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
  @apply grid grid-cols-7 gap-3 px-3 py-2 pb-5;
  p {
    @apply text-xs text-default font-semibold;
  }
  .lc-tooltip-color {
    @apply col-span-1;
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
  .lc-tooltip-step {
    @apply col-span-1;
    @apply w-7;
  }
  .lc-tooltip-value {
    @apply col-span-1;
    @apply w-10 text-left;
  }
  .lc-tooltip-time {
    @apply col-span-2;
  }

  .lc-tooltip-tag {
    @apply col-span-2 truncate;
  }
}
</style>
