package uploader

import (
	"context"
	"fmt"

	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// Result 是一批记录的三态结果：每条输入记录恰好出现在 Confirmed 或 Remaining 之一。
// Confirmed 含确定性失败；Failures 跨两类记录，供诊断使用。
type Result struct {
	Confirmed []*recordv1.Record
	Remaining []*recordv1.Record
	Failures  []Failure
}

// Failure 描述一次失败：Message 面向用户，Records 是受影响的记录，Err 可为 nil。
type Failure struct {
	Message string
	Records []*recordv1.Record
	Err     error
}

// Sender 处理某一类型的整批记录。同一 Uploader 串行调用 Sender。
// Upload 不得修改入参或返回入参之外的记录。
type Sender interface {
	Upload(ctx context.Context, records []*recordv1.Record) Result
}

// Registry 按固定顺序保存 Sender 注册信息，不持有队列或上传协程。
type Registry struct {
	senders map[record.Type]Sender
	order   []record.Type
}

// newRegistry 注册固定分发顺序内的 Sender；未提供或 nil 的类型跳过。
func newRegistry(senders map[record.Type]Sender) *Registry {
	r := &Registry{senders: make(map[record.Type]Sender, len(senders))}
	for _, typ := range dispatchOrder {
		sender, ok := senders[typ]
		if !ok || sender == nil {
			continue
		}
		r.order = append(r.order, typ)
		r.senders[typ] = sender
	}
	return r
}

// has 报告某个类型是否注册了 Sender。
func (r *Registry) has(typ record.Type) bool { return r.senders[typ] != nil }

// get 返回类型对应的 Sender。
func (r *Registry) get(typ record.Type) Sender { return r.senders[typ] }

// len 返回已注册的 Sender 数。
func (r *Registry) len() int { return len(r.senders) }

// safeUpload 兜底 sender 的 panic：整批 remaining，等待下一轮重试。
func safeUpload(ctx context.Context, sender Sender, batch []*recordv1.Record) (result Result) {
	defer func() {
		if recovered := recover(); recovered != nil {
			result = Result{
				Remaining: batch,
				Failures: []Failure{{
					Message: fmt.Sprintf("sender panic: %v", recovered),
					Records: batch,
				}},
			}
		}
	}()
	return sender.Upload(ctx, batch)
}

// conserved 校验每条输入记录恰好出现一次，且没有重复或入参之外的记录。
func conserved(batch []*recordv1.Record, result Result) bool {
	if len(result.Confirmed)+len(result.Remaining) != len(batch) {
		return false
	}
	seen := make(map[*recordv1.Record]struct{}, len(batch))
	for _, rec := range result.Confirmed {
		if _, dup := seen[rec]; dup {
			return false
		}
		seen[rec] = struct{}{}
	}
	for _, rec := range result.Remaining {
		if _, dup := seen[rec]; dup {
			return false
		}
		seen[rec] = struct{}{}
	}
	for _, rec := range batch {
		if _, ok := seen[rec]; !ok {
			return false
		}
	}
	return true
}

// remainingInBatchOrder 按原批序过滤 remaining，避免打乱 FIFO。
func remainingInBatchOrder(batch, remaining []*recordv1.Record) []*recordv1.Record {
	if len(remaining) == 0 {
		return nil
	}
	set := make(map[*recordv1.Record]struct{}, len(remaining))
	for _, rec := range remaining {
		set[rec] = struct{}{}
	}
	kept := make([]*recordv1.Record, 0, len(remaining))
	for _, rec := range batch {
		if _, ok := set[rec]; ok {
			kept = append(kept, rec)
		}
	}
	return kept
}
