// 本文件承载各类 Upsert 的公共入口、两阶段提交与错误映射。
package manager

import (
	"context"
	"errors"
	"fmt"

	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/manager/planner"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

// preparedBatch 封装单批待入队的记录，以及通道预留成功后才执行的会话状态更新（如列定义、LIVE 监听注册）。
type preparedBatch struct {
	records []*recordv1.Record
	apply   func() error
}

// UpsertColumns 定义一批列；已定义的 key 跳过，空列表按 no-op 成功返回。
func (m *Manager) UpsertColumns(ctx context.Context, sessionID string, columns []*columnv1.ColumnRecord) error {
	s, err := m.session(sessionID)
	if err != nil {
		return err
	}
	return s.upsertColumns(ctx, columns)
}

// UpsertScalars 提交一批标量；新 key 自动定义列，重复/过期 step 逐条跳过。
func (m *Manager) UpsertScalars(ctx context.Context, sessionID string, records []*datav1.ScalarRecord) error {
	s, err := m.session(sessionID)
	if err != nil {
		return err
	}
	return s.upsertScalars(ctx, records)
}

// UpsertMedia 提交一批媒体指标；语义与 UpsertScalars 一致。
func (m *Manager) UpsertMedia(ctx context.Context, sessionID string, records []*datav1.MediaRecord) error {
	s, err := m.session(sessionID)
	if err != nil {
		return err
	}
	return s.upsertMedia(ctx, records)
}

// UpsertHistograms 提交一批直方图指标；语义与 UpsertScalars 一致。
func (m *Manager) UpsertHistograms(ctx context.Context, sessionID string, records []*datav1.HistogramRecord) error {
	s, err := m.session(sessionID)
	if err != nil {
		return err
	}
	return s.upsertHistograms(ctx, records)
}

// UpsertLogs 提交一批终端日志，epoch 从 Resume 恢复的行号继续递增。
func (m *Manager) UpsertLogs(ctx context.Context, sessionID string, records []*terminalv1.LogRecord) error {
	s, err := m.session(sessionID)
	if err != nil {
		return err
	}
	return s.upsertLogs(ctx, records)
}

// UpsertSaves 提交一批保存记录：内部类型只做元数据上传，CUSTOM 额外建立
// 本地镜像并按 policy 立即上传或延迟到 Close。
func (m *Manager) UpsertSaves(ctx context.Context, sessionID string, records []*savev1.SaveRecord) error {
	s, err := m.session(sessionID)
	if err != nil {
		return err
	}
	return s.upsertSaves(ctx, records)
}

// commitPrepared 在会话锁内按「状态校验 → 构建记录 → 预留容量 → 更新会话状态 → 提交入队」的顺序提交；
// 预留失败时不会执行 apply，避免污染指标定义或文件监听状态。
func (s *Session) commitPrepared(ctx context.Context, build func() (preparedBatch, error)) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	if RunState(s.state.Load()) != RunStateRunning {
		return ErrNotRunning
	}
	batch, err := build()
	if err != nil {
		return err
	}
	if len(batch.records) == 0 {
		return nil
	}
	ticket, err := s.intake.Reserve(ctx, planner.DemandFor(batch.records, s.cfg))
	if err != nil {
		return upsertError(err)
	}
	if batch.apply != nil {
		if err = batch.apply(); err != nil {
			s.intake.Release(ticket)
			return fmt.Errorf("%w: %v", ErrInvalidRecord, err)
		}
	}
	if err = s.intake.Commit(ticket, batch.records); err != nil {
		return upsertError(err)
	}
	return nil
}

// upsertError 把 intake 的容量错误映射为 manager 的领域错误。
func upsertError(err error) error {
	switch {
	case errors.Is(err, buffer.ErrStoreDegraded):
		return fmt.Errorf("%w: %v", ErrStoreUnavailable, err)
	case errors.Is(err, buffer.ErrOverloaded), errors.Is(err, buffer.ErrBatchTooLarge):
		return fmt.Errorf("%w: %v", ErrOverloaded, err)
	case errors.Is(err, buffer.ErrClosed):
		return fmt.Errorf("%w: %w", ErrNotRunning, err)
	default:
		return err
	}
}

// commitChunks 按单批上限切分请求并依次提交，适用于支持按 key/step 幂等去重的记录。
func commitChunks[T any](ctx context.Context, s *Session, items []T, prepare func([]T) (preparedBatch, error)) error {
	size := s.chunkSize()
	for start := 0; start < len(items); start += size {
		chunk := items[start:min(start+size, len(items))]
		if err := s.commitPrepared(ctx, func() (preparedBatch, error) {
			return prepare(chunk)
		}); err != nil {
			return err
		}
	}
	return nil
}

func (s *Session) chunkSize() int {
	if s.cfg.RecordBatch <= 0 {
		return 1
	}
	return int(s.cfg.RecordBatch)
}

// checkBatchLimit 校验单批记录数不超过上限。日志与保存记录缺乏去重键，
// 若拆批后部分成功会导致重试时产生重复数据，因此超限时整批拒绝。
func (s *Session) checkBatchLimit(n int) error {
	limit := s.chunkSize()
	if n > limit {
		return fmt.Errorf("%w: %d records exceed record_batch %d", ErrRequestTooLarge, n, limit)
	}
	return nil
}
