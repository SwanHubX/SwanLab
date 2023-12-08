import moment from 'moment'
import { t } from '@swanlab-vue/i18n'
/**
 * 时间格式转换，转换为相对时间，转换结果由i18n提供
 * @param { number } timestamp - 时间戳 一个时间的格林威治时间数值
 * @returns { string } "5 months ago"
 */
export function transTime(timestamp) {
  moment.updateLocale('zh', {
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
  return moment(timestamp).fromNow()
}
