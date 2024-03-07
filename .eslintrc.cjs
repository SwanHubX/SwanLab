module.exports = {
  env: {
    browser: true,
    es2021: true,
    node: true
  },
  extends: [
    'plugin:vue/base',
    'eslint:recommended',
    'plugin:vue/vue3-essential',
    './.eslintrc-auto-import.json'
    // "plugin:vue/vue3-strongly-recommended"
  ],
  overrides: [],
  parserOptions: {
    ecmaVersion: 'latest',
    sourceType: 'module'
  },
  plugins: ['vue'],
  rules: {}
}
