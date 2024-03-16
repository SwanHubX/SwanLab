/**
 * 在此处封装一些图表可共用的函数逻辑
 */

/**
 * 将refrence字段映射为x轴字段和显示的x轴名称
 */
export const refrence2XField = {
  step: { xField: 'index', xTitle: 'Step' }
}

export const transparentColor = (color, opacity = 0.1) => {
  // 将十六进制增加透明度通道
  opacity = Math.floor(opacity * 255)
  // opacity转为十六进制
  const opacityHex = opacity.toString(16)
  return color + opacityHex
}
