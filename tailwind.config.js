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
        full: '9999'
      },
      // 背景颜色
      backgroundColor: {
        default: 'var(--background-default, #fff)',
        dimmer: 'var(--background-dimmer, #f6f8fA)'
      },
      // 边框颜色
      borderColor: {
        default: 'var(--outline-default, #d1d7dd)'
      },
      // 文字颜色
      textColor: {
        default: 'var(--foreground-default, #000)',
        dimmer: 'var(--foreground-dimmer, #555)',
        dimmest: 'var(--foreground-dimmest, #b5b5b5)'
      },
      // 基础颜色
      colors: {
        // 危险颜色
        'negative-dimmest': 'var(--nagative-dimmest, #eddcdc)',
        'negative-default': 'var(--negative-default, #bb4c51)',
        // 成功颜色
        'positive-dimmer': 'var(--positive-dimmer, #37733c)',
        'positive-default': 'var(--positive-default, #4e9553)',
        'positive-higher': 'var(--positive-higher, #cff1d8)',
        'positive-highest': 'var(--positive-highest, #d5e4d9)',
        // 标准颜色
        'primary-dimmest': 'var(--primary-dimmest, #d6e8ff)',
        'primary-default': 'var(--primary-default, #1c74dd)'
      }
    }
  },
  plugins: []
}
