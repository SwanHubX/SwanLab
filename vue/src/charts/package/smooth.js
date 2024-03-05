/**
 * 数据平滑处理
 */

/**
 * 数据平滑处理，传入数据和平滑方法
 * @param { Array } data 数据
 * @param { Object } method 平滑方法
 * @param { Number } method.id 平滑方法id
 * @param { Number} method.value 平滑方法参数
 * @returns { Array } 平滑后的数据
 */
export default function smooth(data, method) {
  console.log('smooth', data, method)
  data = data.map((item) => {
    return {
      ...item,
      data: item.data * 0.5
    }
  })
  return data
}
