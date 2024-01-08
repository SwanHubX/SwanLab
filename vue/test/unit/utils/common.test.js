import { uuid } from '@swanlab-vue/utils/common'
import { describe, test, expect } from 'vitest'

describe('common.js utils', () => {
  test('generate uuid', () => {
    expect(uuid()).toHaveLength(10)
  })
})
