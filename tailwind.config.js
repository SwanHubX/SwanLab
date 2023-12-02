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
        default: 'var(--background-default)',
        dimmer: 'var(--background-dimmer)'
      },
      // 边框颜色
      borderColor: {
        default: 'var(--outline-default)'
      },
      // 文字颜色
      textColor: {
        default: 'var(--foreground-default)',
        dimmer: 'var(--foreground-dimmer)',
        dimmest: 'var(--foreground-dimmest)',
      },
      // 基础颜色
      colors: {
        // 危险颜色
        'negative-dimmest': 'var(--nagative-dimmest)',
        'negative-default': 'var(--negative-default)',
        // 成功颜色
        'positive-dimmer': 'var(--positive-dimmer)',
        'positive-default': 'var(--positive-default)',
        'positive-higher': 'var(--positive-higher)',
        'positive-highest': 'var(--positive-highest)',
        // 标准颜色
        'primary-dimmest': 'var(--primary-dimmest)',
        'primary-default': 'var(--primary-default)',
      }
    }
  },
  plugins: []
}
