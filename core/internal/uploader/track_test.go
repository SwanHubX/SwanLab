package uploader

import (
	"context"
	"sync"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// senderFunc 是 Sender 的函数适配器。
type senderFunc func(context.Context, []*recordv1.Record) Result

func (f senderFunc) Upload(ctx context.Context, records []*recordv1.Record) Result {
	return f(ctx, records)
}

// scriptedSender 按调用序号返回预置结果，超出脚本后重复最后一条。
type scriptedSender struct {
	mu      sync.Mutex
	scripts []func([]*recordv1.Record) Result
	calls   [][]*recordv1.Record
}

func (s *scriptedSender) Upload(_ context.Context, records []*recordv1.Record) Result {
	s.mu.Lock()
	index := len(s.calls)
	s.calls = append(s.calls, records)
	var step func([]*recordv1.Record) Result
	if index < len(s.scripts) {
		step = s.scripts[index]
	} else if len(s.scripts) > 0 {
		step = s.scripts[len(s.scripts)-1]
	}
	s.mu.Unlock()

	if step == nil {
		return Result{Confirmed: records}
	}
	return step(records)
}

// quietThrottle 返回不输出任何内容的告警节流器。
func quietThrottle() *WarningThrottle {
	throttle := newWarningThrottle(time.Now, time.Hour)
	throttle.warnf = func(string) {}
	throttle.infof = func(string) {}
	return throttle
}

// retryNoJitter 返回无抖动的退避：每次失败固定等待 base。
func retryNoJitter(base time.Duration) RetryPolicy {
	return RetryPolicy{Base: base, Multiplier: 2, Max: base, JitterLow: 1, JitterHigh: 1}
}

func waitSignal(t *testing.T, ch <-chan struct{}, what string) {
	t.Helper()
	select {
	case <-ch:
	case <-time.After(2 * time.Second):
		t.Fatalf("timeout waiting for %s", what)
	}
}

func TestRegistryKeepsDispatchOrder(t *testing.T) {
	senders := map[record.Type]Sender{}
	for _, typ := range dispatchOrder {
		senders[typ] = senderFunc(successUpload)
	}
	senders[record.TypeStart] = senderFunc(successUpload)
	registry := newRegistry(senders)

	if registry.len() != len(dispatchOrder) || len(registry.order) != len(dispatchOrder) {
		t.Fatalf("registered = %d/%d, want %d", registry.len(), len(registry.order), len(dispatchOrder))
	}
	for i, typ := range dispatchOrder {
		if registry.order[i] != typ || !registry.has(typ) || registry.get(typ) == nil {
			t.Fatalf("registration[%d] = %v, want %s", i, registry.order[i], typ)
		}
	}
	if registry.has(record.TypeStart) || registry.get(record.TypeStart) != nil {
		t.Fatal("type outside dispatchOrder must not be registered")
	}
}

func TestRegistrySkipsNilSender(t *testing.T) {
	registry := newRegistry(map[record.Type]Sender{record.TypeScalar: nil})
	if registry.len() != 0 || registry.has(record.TypeScalar) {
		t.Fatal("nil sender must be skipped")
	}
}

func TestSafeUploadRecoversSenderPanic(t *testing.T) {
	records := []*recordv1.Record{newScalar(1, "a"), newScalar(2, "b")}
	result := safeUpload(context.Background(), senderFunc(func(context.Context, []*recordv1.Record) Result {
		panic("boom")
	}), records)
	if len(result.Remaining) != 2 || len(result.Confirmed) != 0 {
		t.Fatalf("result = %+v, want whole batch remaining", result)
	}
	if len(result.Failures) != 1 || result.Failures[0].Message != "sender panic: boom" {
		t.Fatalf("failures = %+v, want panic failure", result.Failures)
	}
}

func TestRemainingKeepsBatchOrder(t *testing.T) {
	a, b, c := newScalar(1, "a"), newScalar(2, "b"), newScalar(3, "c")
	kept := remainingInBatchOrder(
		[]*recordv1.Record{a, b, c},
		[]*recordv1.Record{c, a},
	)
	if len(kept) != 2 || kept[0] != a || kept[1] != c {
		t.Fatalf("remaining order = %v, want [a c]", kept)
	}
}

func TestConserved(t *testing.T) {
	a, b := newScalar(1, "a"), newScalar(2, "b")
	batch := []*recordv1.Record{a, b}

	cases := []struct {
		name   string
		result Result
		want   bool
	}{
		{"all confirmed", Result{Confirmed: batch}, true},
		{"all remaining", Result{Remaining: batch}, true},
		{"split", Result{Confirmed: []*recordv1.Record{a}, Remaining: []*recordv1.Record{b}}, true},
		{"dropped", Result{Confirmed: []*recordv1.Record{a}}, false},
		{"duplicate", Result{Confirmed: batch, Remaining: []*recordv1.Record{a}}, false},
		{"extra record", Result{Confirmed: append([]*recordv1.Record{a, b}, newScalar(3, "c"))}, false},
		{"nothing", Result{}, false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := conserved(batch, tc.result); got != tc.want {
				t.Fatalf("conserved = %v, want %v", got, tc.want)
			}
		})
	}
}
