// Package uploadertest 提供 uploader.Sender 的测试替身。
package uploadertest

import (
	"context"
	"sync"

	"github.com/swanhubx/swanlab/core/internal/uploader"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// Sender 是可编排的假 sender：记录调用批次，并按序号返回预置结果（并发安全）。
type Sender struct {
	mu    sync.Mutex
	steps []func([]*recordv1.Record) uploader.Result
	calls [][]*recordv1.Record
}

// New 构造全部确认的假 sender。
func New() *Sender { return &Sender{} }

// AlwaysRetry 返回永远可重试的 sender。
func AlwaysRetry() *Sender {
	return New().On(func(records []*recordv1.Record) uploader.Result {
		return uploader.Result{Remaining: records}
	})
}

// On 追加脚本步骤；超出后重复最后一步，无脚本时全部确认。
func (s *Sender) On(step func([]*recordv1.Record) uploader.Result) *Sender {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.steps = append(s.steps, step)
	return s
}

// Upload 实现 uploader.Sender。
func (s *Sender) Upload(_ context.Context, records []*recordv1.Record) uploader.Result {
	s.mu.Lock()
	index := len(s.calls)
	s.calls = append(s.calls, append([]*recordv1.Record(nil), records...))
	var step func([]*recordv1.Record) uploader.Result
	if index < len(s.steps) {
		step = s.steps[index]
	} else if len(s.steps) > 0 {
		step = s.steps[len(s.steps)-1]
	}
	s.mu.Unlock()

	if step == nil {
		return uploader.Result{Confirmed: records}
	}
	return step(records)
}

// Calls 返回调用批次快照。
func (s *Sender) Calls() [][]*recordv1.Record {
	s.mu.Lock()
	defer s.mu.Unlock()
	return append([][]*recordv1.Record(nil), s.calls...)
}

// Count 返回调用次数。
func (s *Sender) Count() int {
	s.mu.Lock()
	defer s.mu.Unlock()
	return len(s.calls)
}
