// Package budget 提供上传管线的会话级容量预算（条数与字节双限）。
//
// 预算覆盖记录从入队到上传终结的整个生命周期：记录在组件间转移时不释放
// 配额，只有同号覆盖或确认终结才归还。上限 <= 0 表示该维度不限制。
package budget

import (
	"context"
	"errors"
	"math"
	"sync"
)

// ErrTooLarge 表示单次申请量本身超过上限，等待也不会成功。
var ErrTooLarge = errors.New("budget: request exceeds the configured limit")

// Limits 是条数与字节上限，<=0 表示该维度不限制。
type Limits struct {
	Records int64
	Bytes   int64
}

// Unlimited 报告两个维度是否都不限制。
func (l Limits) Unlimited() bool { return l.Records <= 0 && l.Bytes <= 0 }

// Budget 是并发安全的容量账本。零值不可用，应由 New 构造。
type Budget struct {
	mu          sync.Mutex
	limits      Limits
	usedRecords int64
	usedBytes   int64
	// wait 在每次归还后关闭并重建，用于唤醒阻塞中的 Acquire。
	wait chan struct{}
}

// New 创建容量账本；limits 的负值按不限制处理。
func New(limits Limits) *Budget {
	if limits.Records < 0 {
		limits.Records = 0
	}
	if limits.Bytes < 0 {
		limits.Bytes = 0
	}
	return &Budget{limits: limits, wait: make(chan struct{})}
}

// Limits 返回构造时的上限。
func (b *Budget) Limits() Limits { return b.limits }

// Usage 返回当前占用的条数与字节数。
func (b *Budget) Usage() (records, bytes int64) {
	b.mu.Lock()
	defer b.mu.Unlock()
	return b.usedRecords, b.usedBytes
}

// TryAcquire 尝试立即占用配额；额度不足或单次申请超过上限时返回 false，
// 且不改变任何状态。
func (b *Budget) TryAcquire(records, bytes int64) bool {
	checkAmount(records, bytes)
	b.mu.Lock()
	defer b.mu.Unlock()
	if b.tooLarge(records, bytes) || !b.fitsLocked(records, bytes) {
		return false
	}
	b.usedRecords += records
	b.usedBytes += bytes
	return true
}

// Acquire 占用配额，额度不足时阻塞到有归还或 ctx 结束。
// 单次申请量超过上限时立即返回 ErrTooLarge；不保证 FIFO。
func (b *Budget) Acquire(ctx context.Context, records, bytes int64) error {
	checkAmount(records, bytes)
	if b.tooLarge(records, bytes) {
		return ErrTooLarge
	}
	for {
		b.mu.Lock()
		if b.fitsLocked(records, bytes) {
			b.usedRecords += records
			b.usedBytes += bytes
			b.mu.Unlock()
			return nil
		}
		wait := b.wait
		b.mu.Unlock()

		select {
		case <-ctx.Done():
			return ctx.Err()
		case <-wait:
		}
	}
}

// Release 归还配额并唤醒等待者。
// 超量归还属于调用约定错误，会 panic，避免账本静默失衡。
func (b *Budget) Release(records, bytes int64) {
	checkAmount(records, bytes)
	b.mu.Lock()
	defer b.mu.Unlock()
	if records > b.usedRecords || bytes > b.usedBytes {
		panic("budget: released more than reserved")
	}
	b.usedRecords -= records
	b.usedBytes -= bytes
	if records > 0 || bytes > 0 {
		b.wakeLocked()
	}
}

// Adjust 按实际用量修正一次预留：正数补扣、负数退还。
// 补扣可能使占用短暂超过上限；退款超量会 panic。
func (b *Budget) Adjust(records, bytes int64) {
	b.mu.Lock()
	defer b.mu.Unlock()
	if records < 0 && -records > b.usedRecords || bytes < 0 && -bytes > b.usedBytes {
		panic("budget: adjusted more than reserved")
	}
	b.usedRecords += records
	b.usedBytes += bytes
	if records < 0 || bytes < 0 {
		b.wakeLocked()
	}
}

// Full 报告任一有限维度是否已经用满。
func (b *Budget) Full() bool {
	b.mu.Lock()
	defer b.mu.Unlock()
	return b.limits.Records > 0 && b.usedRecords >= b.limits.Records ||
		b.limits.Bytes > 0 && b.usedBytes >= b.limits.Bytes
}

// Ratio 返回占用率：条数与字节取较大值，无限维度不参与计算。
// 两个维度都不限制时返回 0；经 Adjust 后可能短暂大于 1。
func (b *Budget) Ratio() float64 {
	b.mu.Lock()
	defer b.mu.Unlock()
	var ratio float64
	if b.limits.Records > 0 {
		ratio = math.Max(ratio, float64(b.usedRecords)/float64(b.limits.Records))
	}
	if b.limits.Bytes > 0 {
		ratio = math.Max(ratio, float64(b.usedBytes)/float64(b.limits.Bytes))
	}
	return ratio
}

// wakeLocked 唤醒所有等待者；调用方需持有 mu。
func (b *Budget) wakeLocked() {
	close(b.wait)
	b.wait = make(chan struct{})
}

// fitsLocked 判断当前占用加上本次申请是否仍在上限内；调用方需持有 mu。
func (b *Budget) fitsLocked(records, bytes int64) bool {
	if b.limits.Records > 0 && b.usedRecords+records > b.limits.Records {
		return false
	}
	if b.limits.Bytes > 0 && b.usedBytes+bytes > b.limits.Bytes {
		return false
	}
	return true
}

// tooLarge 判断单次申请量是否已超过上限；limits 不可变，无需加锁。
func (b *Budget) tooLarge(records, bytes int64) bool {
	return b.limits.Records > 0 && records > b.limits.Records ||
		b.limits.Bytes > 0 && bytes > b.limits.Bytes
}

func checkAmount(records, bytes int64) {
	if records < 0 || bytes < 0 {
		panic("budget: amount must be non-negative")
	}
}
