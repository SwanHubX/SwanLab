// Package counter 提供并发安全的自增计数器。
//
// Record 编号（builder）与控制台行号（epoch）都由它分配；Inc/Dec 返回
// 执行后的新值，可直接作为 Record.num 使用。
package counter

import "sync/atomic"

// Counter 是并发安全的整数计数器。
//
// 零值可用，初始值为 0。
type Counter struct {
	v atomic.Int64
}

// New 创建初始值为 initial 的计数器。
func New(initial int64) *Counter {
	c := &Counter{}
	c.v.Store(initial)
	return c
}

// Inc 计数加 delta（默认路径 delta=1），返回增加后的值。
func (c *Counter) Inc() int64 { return c.v.Add(1) }

// IncBy 计数加 delta，返回增加后的值；delta 可为负数。
func (c *Counter) IncBy(delta int64) int64 { return c.v.Add(delta) }

// Dec 计数减 1，返回减少后的值。
func (c *Counter) Dec() int64 { return c.v.Add(-1) }

// DecBy 计数减 delta，返回减少后的值；delta 可为负数。
func (c *Counter) DecBy(delta int64) int64 { return c.v.Add(-delta) }

// Value 返回当前计数值。
func (c *Counter) Value() int64 { return c.v.Load() }

// Reset 将计数器重置为 value。
func (c *Counter) Reset(value int64) { c.v.Store(value) }

// GetAndReset 返回当前值并将计数器重置为 value。
func (c *Counter) GetAndReset(value int64) int64 { return c.v.Swap(value) }
