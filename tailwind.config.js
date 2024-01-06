const length = {
  0.25: '1px',
  4.5: '18px',
  6.5: '26px',
  7.5: '30px',
  15: '60px',
  66: '264px'
}

/** @type {import('tailwindcss').Config} */
export default {
  content: ['./index.html', './vue/src/**/*.{vue,js,ts,jsx,tsx}'],
  theme: {
    extend: {
      margin: {
        ...length
      },
      padding: {
        ...length
      },
      width: {
        ...length
      },
      height: {
        ...length
      },
      minHeight: {
        ...length
      },
      zIndex: {
        // 最顶层，给message用
        message: '9999',
        // 其余最顶层
        full: '9998'
      },
      // 背景颜色
      backgroundColor: {
        overlay: `var(--background-overlay, rgba(15, 21, 36, 0.8))`,
        default: 'var(--background-default, #fff)',
        higher: 'var(--background-higher, #f6f8fA)',
        highest: 'var(--background-highest, #ebebec)'
      },
      // 边框颜色
      borderColor: {
        default: 'var(--outline-default, #d1d7dd)'
      },
      // 文字颜色
      textColor: {
        default: 'var(--foreground-default, #000)',
        dimmer: 'var(--foreground-dimmer, #555)',
        dimmest: 'var(--foreground-dimmest, #b5b5b5)',
        // 复杂文字高亮
        'white-default': 'var(--accent-white-default, #fff)'
      }
    },
    // 基础颜色
    colors: {
      transparent: 'transparent',
      current: 'currentColor',
      /* 成功颜色 */
      'positive-highest': 'var(--positive-highest, #004d0d)',
      'positive-higher': 'var(--positive-higher, #036e15)',
      'positive-default': 'var(--positive-default, #00a11b)',
      'positive-dimmer': 'var(--positive-dimmer, #3cc954)',
      'positive-dimmest': 'var(--positive-dimmest, #97e896)',
      /* 标准颜色 */
      'primary-highest': 'var(--primary-highest, #004182)',
      'primary-higher': 'var(--primary-higher, #005cb8)',
      'primary-default': 'var(--primary-default, #0f87ff)',
      'primary-dimmer': 'var(--primary-dimmer, #74b9ff)',
      'primary-dimmest': 'var(--primary-dimmest, #b2d9ff)',
      /* 警告颜色 */
      'warning-dimmer': 'var(--warning-dimmer, #ff9933)',
      'warning-dimmest': 'var(--warning-dimmest, #ffad5a)',
      /* 危险颜色 */
      'negative-highest': 'var(--negative-highest, #8a0000)',
      'negative-higher': 'var(--negative-higher, #c20a0a)',
      'negative-default': 'var(--negative-default, #fa4b4b)',
      'negative-dimmest': 'var(--negative-dimmest, #ff9494)',
      'negative-dimmer': 'var(--negative-dimmer, #ffc7c7)'
    }
  },
  plugins: []
}
