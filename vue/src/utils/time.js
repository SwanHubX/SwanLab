import moment from 'moment'
import { t } from '@swanlab-vue/i18n'
/**
 * 时间格式转换，转换为相对时间，转换结果由i18n提供
 * 相对时间只管显示，不管时区，如果需要时区转换，请使用convertUtcToLocal
 * @param { number } timestamp - 时间戳,单位为秒
 * @returns { string } "5 months ago"
 */
export function transTime(time) {
  moment.updateLocale('local', {
    relativeTime: {
      future: t('moment.relativeTime.future'), // 未来的算在之前
      past: t('moment.relativeTime.past'),
      s: t('moment.relativeTime.s'),
      ss: t('moment.relativeTime.ss'),
      m: t('moment.relativeTime.m'),
      mm: t('moment.relativeTime.mm'),
      h: t('moment.relativeTime.h'),
      hh: t('moment.relativeTime.hh'),
      d: t('moment.relativeTime.d'),
      dd: t('moment.relativeTime.dd'),
      w: t('moment.relativeTime.w'),
      ww: t('moment.relativeTime.ww'),
      M: t('moment.relativeTime.M'),
      MM: t('moment.relativeTime.MM'),
      y: t('moment.relativeTime.y'),
      yy: t('moment.relativeTime.yy')
    }
  })
  // 根据传入的时间和时区，转换为当前浏览器下的时间

  return moment(time).fromNow()
}

/**
 * 将 UTC 时间转换为当前浏览器时区的时间。
 *
 * @param {string} utcTime - UTC 时间字符串，使用 ISO 8601 格式。
 * @returns {Date} 格式化的本地时间
 * @throws {Error} 如果输入的时间字符串格式不正确。
 *
 * @example
 * // 使用示例
 * const utcTime = "2023-01-01T12:00:00Z";
 * const localTime = convertUtcToLocal(utcTime);
 * console.log("本地时间:", localTime);
 */
export function convertUtcToLocal(utcTime) {
  // 创建 UTC 时间的 Date 对象
  var utcDate = new Date(utcTime)

  // 获取本地时区的时间偏移（分钟）
  var localTimeZoneOffset = new Date().getTimezoneOffset()

  // 计算 UTC 时间与本地时间的时间差（毫秒）
  var timeZoneDifference = localTimeZoneOffset * 60 * 1000

  // 计算本地时间
  var localTime = new Date(utcDate.getTime() - timeZoneDifference)

  // 检查结果是否有效
  if (isNaN(localTime.getTime())) {
    throw new Error('Invalid UTC time format')
  }

  return localTime
}

/**
 * 格式化时间，接收一个时间戳，转化当前浏览器当前时区下的时间，输入时间为UTC时间
 * 例如：上海时区下，2023-12-04T04:44:33.026550Z转换为"2023/12/04 12:44:33"
 * @param {string} time - 时间字符串
 */

export const formatTime = (time) => {
  // 判断是否携带时区信息
  const zone = typeof time === 'string' && time.endsWith('Z')
  let { year, month, day, hour, minute, second } = getTimes(time, zone)
  // 如果minute是个位数，转换成两位数
  if (month < 10) month = '0' + month

  if (day < 10) day = '0' + day

  if (hour < 10) hour = '0' + hour

  if (minute < 10) minute = '0' + minute

  if (second < 10) second = '0' + second

  return `${year}/${month}/${day} ${hour}:${minute}:${second}`
}

/**
 * 获取年月日时分秒
 * @param {string} time - 时间字符串
 * @param {boolean} zone - 是否需要时区转换，为true意味着携带时区信息，不需要时区转换
 */
export const getTimes = (time, zone = true) => {
  let localDate
  if (!zone) {
    const date = new Date(time)
    const timezoneOffset = new Date().getTimezoneOffset()
    const localTime = date.getTime() - timezoneOffset * 60 * 1000
    localDate = new Date(localTime)
  } else {
    localDate = new Date(time)
  }
  const year = localDate.getFullYear()
  let month = localDate.getMonth() + 1
  let day = localDate.getDate()
  let hour = localDate.getHours()
  let minute = localDate.getMinutes()
  let second = localDate.getSeconds()

  return {
    year,
    month,
    day,
    hour,
    minute,
    second
  }
}

/**
 * 获取实验的持续时间
 * @param {experiment} experiment 实验对象
 * @returns {string} 实验持续时间
 */
export const getDuration = (experiment) => {
  if (!experiment) return null
  const time1 = new Date(experiment.create_time)
  const currentTime = new Date()
  const time2 = experiment.status === 0 ? new Date(currentTime.getTime()) : new Date(experiment.finish_time)

  if (isNaN(time1.getTime()) || isNaN(time2.getTime())) {
    // 处理无效日期的情况
    return 'Invalid date'
  }

  const timeDifference = Math.abs(time2 - time1)

  const seconds = Math.floor(timeDifference / 1000) % 60
  const minutes = Math.floor(timeDifference / (1000 * 60)) % 60
  const hours = Math.floor(timeDifference / (1000 * 60 * 60)) % 24
  const days = Math.floor(timeDifference / (1000 * 60 * 60 * 24))

  const formattedTime = []

  if (days > 0) {
    formattedTime.push(`${days}${t('experiment.index.header.experiment_infos.time.day')}`)
  }

  if (hours > 0) {
    formattedTime.push(`${hours}${t('experiment.index.header.experiment_infos.time.hour')}`)
  }

  if (minutes > 0) {
    formattedTime.push(`${minutes}${t('experiment.index.header.experiment_infos.time.minute')}`)
  }

  if (seconds > 0) {
    formattedTime.push(`${seconds}${t('experiment.index.header.experiment_infos.time.second')}`)
  }

  return formattedTime.join('') || 'less than 1s'
}
