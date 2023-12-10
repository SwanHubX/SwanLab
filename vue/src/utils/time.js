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
