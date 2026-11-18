package syncer

import (
	"errors"
	"fmt"
	"io"

	"google.golang.org/protobuf/proto"
	"google.golang.org/protobuf/types/known/timestamppb"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/store"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

func (s *syncSession) scanLoop() {
	defer close(s.scanDone)

	reader, err := store.Open(s.runFile)
	if err != nil {
		s.finishScan(fmt.Errorf("failed to open run file %s: %w", s.runFile, err))
		return
	}
	defer func() { _ = reader.Close() }()

	// 首条 StartRecord 已在 DeliverSyncStart 校验过；这里重验一次，避免源文件被同名替换后把新数据写进旧目标。
	if verr := s.verifyStartRecord(reader); verr != nil {
		s.finishScan(verr)
		return
	}

	for {
		if s.ctx.Err() != nil {
			return
		}
		payload, scanErr := reader.Scan()
		if scanErr != nil {
			switch {
			case errors.Is(scanErr, io.EOF):
			case errors.Is(scanErr, store.ErrTruncated), errors.Is(scanErr, store.ErrCorrupt):
				console.Warningf("sync: stopped scanning %s due to store corruption: %v", s.runFile, scanErr)
			default:
				// 真实 I/O 失败不能当成“正常扫描到尾”，必须传给 Confirm 拒绝成功。
				s.finishScan(fmt.Errorf("failed to scan run file %s: %w", s.runFile, scanErr))
				return
			}
			break
		}

		rec := &recordv1.Record{}
		if uerr := proto.Unmarshal(payload, rec); uerr != nil {
			console.Warningf("sync: failed to parse record payload, skipping: %v", uerr)
			continue
		}
		if !s.handleRecord(rec) {
			return
		}
	}

	s.finishScan(nil)
}

// verifyStartRecord 消费首条 StartRecord，并校验它与 Start 阶段读到的 run id 一致。
func (s *syncSession) verifyStartRecord(reader *store.Reader) error {
	payload, err := reader.Scan()
	if err != nil {
		if errors.Is(err, io.EOF) {
			// 这些常量会直接作为 Confirm 的业务消息返回给用户，故意保留首字母大写。
			return fmt.Errorf(missingEmptyRunFileMsg, s.runFile) //nolint:staticcheck // ST1005
		}
		return fmt.Errorf("failed to read run file %s: %w", s.runFile, err)
	}
	rec := &recordv1.Record{}
	if err := proto.Unmarshal(payload, rec); err != nil {
		return fmt.Errorf("failed to parse run file %s: %w", s.runFile, err)
	}
	start := rec.GetStart()
	if start == nil {
		return errors.New(missingStartRecordMsg) //nolint:staticcheck // ST1005
	}
	if want := s.startRecord.GetId(); want != "" && start.GetId() != want {
		return fmt.Errorf("run file %s was replaced: start id %q does not match %q", s.runFile, start.GetId(), want)
	}
	return nil
}

func (s *syncSession) finishScan(scanErr error) {
	if s.ctx.Err() != nil {
		return
	}
	s.mu.Lock()
	s.scanErr = scanErr
	s.mu.Unlock()
	if scanErr == nil {
		s.emitCrashLog()
	}
	close(s.inCh)
	s.uploader.RequestFinish()
}

func (s *syncSession) emitCrashLog() {
	if s.filter.Enabled() || s.finishRecord == nil {
		return
	}
	if s.finishRecord.GetState() == runv1.RunState_RUN_STATE_FINISHED {
		return
	}
	line := firstNonEmpty(s.finishRecord.GetError(), finishUnknownError)
	rec := &recordv1.Record{
		Num:       s.negCounter.Dec(),
		Timestamp: timestamppb.Now(),
		RecordType: &recordv1.Record_Log{
			Log: &terminalv1.LogRecord{
				Line:      line,
				Level:     terminalv1.LogLevel_LOG_LEVEL_ERROR,
				Timestamp: s.finishRecord.GetFinishedAt(),
				Epoch:     s.epoch.Inc(),
			},
		},
	}
	_ = s.enqueue(rec)
}

func (s *syncSession) handleRecord(rec *recordv1.Record) bool {
	switch record.TypeOf(rec) {
	case record.TypeStart:
		return true
	case record.TypeFinish:
		s.handleFinish(rec.GetFinish())
		return true
	case record.TypeLog:
		return s.handleLog(rec)
	case record.TypeColumn:
		return s.handleColumn(rec)
	case record.TypeScalar:
		return s.handleScalar(rec)
	case record.TypeMedia:
		return s.handleMedia(rec)
	case record.TypeHistogram:
		return s.handleHistogram(rec)
	case record.TypeSave:
		return s.enqueue(rec)
	default:
		return true
	}
}

func (s *syncSession) handleFinish(finish *runv1.FinishRecord) {
	if finish == nil {
		return
	}
	if s.finishRecord == nil {
		s.finishRecord = finish
		return
	}
	console.Warning("Multiple finish records were found in the run file. SwanLab will use the first one.")
}

func (s *syncSession) handleLog(rec *recordv1.Record) bool {
	logRec := rec.GetLog()
	if logRec == nil {
		return true
	}
	if s.filter.Enabled() {
		if !s.filter.HasTime() {
			return true
		}
		if !s.filter.MatchTime(PickTimestamp(logRec.GetTimestamp(), rec.GetTimestamp())) {
			return true
		}
		if logRec.GetEpoch() > s.epoch.Value() {
			s.epoch.Reset(logRec.GetEpoch())
		}
		return s.enqueue(rec)
	}
	if logRec.GetEpoch() > s.epoch.Value() {
		s.epoch.Reset(logRec.GetEpoch())
		return s.enqueue(rec)
	}
	console.Debugf("Log at epoch %d was skipped syncing because it is too old.", logRec.GetEpoch())
	return true
}

func (s *syncSession) handleColumn(rec *recordv1.Record) bool {
	col := rec.GetColumn()
	if col == nil || col.GetColumnKey() == "" {
		return true
	}
	key := col.GetColumnKey()
	if s.metrics.Defined(key) {
		console.Debugf("Column '%s' was skipped syncing because it already exists.", key)
		return true
	}
	if err := s.metrics.Define(key, col, -1, columnKind(col.GetColumnType())); err != nil {
		console.Warningf("sync: define column %q failed: %v", key, err)
		return true
	}
	return s.enqueue(rec)
}

func (s *syncSession) handleScalar(rec *recordv1.Record) bool {
	data := rec.GetScalar()
	if data == nil {
		return true
	}
	return s.handleMetric(
		data.GetKey(),
		data.GetStep(),
		data.GetType(),
		PickTimestamp(data.GetTimestamp(), rec.GetTimestamp()),
		rec,
		func() metrics.Precheck { return s.metrics.PrecheckScalar(data) },
		func(col *columnv1.ColumnRecord) (bool, error) { return s.metrics.CommitScalar(data, col) },
	)
}

func (s *syncSession) handleMedia(rec *recordv1.Record) bool {
	data := rec.GetMedia()
	if data == nil {
		return true
	}
	return s.handleMetric(
		data.GetKey(),
		data.GetStep(),
		data.GetType(),
		PickTimestamp(data.GetTimestamp(), rec.GetTimestamp()),
		rec,
		func() metrics.Precheck { return s.metrics.PrecheckMedia(data) },
		func(col *columnv1.ColumnRecord) (bool, error) { return s.metrics.CommitMedia(data, col) },
	)
}

func (s *syncSession) handleHistogram(rec *recordv1.Record) bool {
	data := rec.GetHistogram()
	if data == nil {
		return true
	}
	return s.handleMetric(
		data.GetKey(),
		data.GetStep(),
		data.GetType(),
		PickTimestamp(data.GetTimestamp(), rec.GetTimestamp()),
		rec,
		func() metrics.Precheck { return s.metrics.PrecheckHistogram(data) },
		func(col *columnv1.ColumnRecord) (bool, error) { return s.metrics.CommitHistogram(data, col) },
	)
}

func (s *syncSession) handleMetric(
	key string,
	step int64,
	typ columnv1.ColumnType,
	ts *timestamppb.Timestamp,
	rec *recordv1.Record,
	precheckFn func() metrics.Precheck,
	commitFn func(*columnv1.ColumnRecord) (bool, error),
) bool {
	if key == "" {
		return true
	}
	if s.filter.Enabled() && !s.filter.MatchMetric(step, ts) {
		return true
	}

	pre := precheckFn()
	var autoCol *columnv1.ColumnRecord
	if pre.NeedsColumn {
		autoCol = builder.BuildAutoColumn(s.cfg.SectionRule, key, typ)
		colRec := &recordv1.Record{
			Num:       s.negCounter.Dec(),
			Timestamp: timestamppb.Now(),
			RecordType: &recordv1.Record_Column{
				Column: autoCol,
			},
		}
		if !s.enqueue(colRec) {
			return false
		}
	}

	accepted, err := commitFn(autoCol)
	if err != nil {
		console.Warningf("sync: metric %q at step %d skipped: %v", key, step, err)
		return true
	}
	if !accepted {
		console.Debugf("Metric '%s' at step %d was skipped syncing because it is duplicate or too old.", key, step)
		return true
	}
	return s.enqueue(rec)
}

func (s *syncSession) enqueue(rec *recordv1.Record) bool {
	// 入队前按字节配额阻塞预留并登记总量，确保流式读盘不堆积大量内联载荷，
	// 单条超限时走 Adjust 放行避免阻塞，且不会出现 uploaded > total 的瞬态倒挂。
	size := record.Size(rec)
	if err := s.budget.Acquire(s.ctx, 1, size); err != nil {
		if !errors.Is(err, budget.ErrTooLarge) {
			return false
		}
		s.budget.Adjust(1, size)
	}
	s.tracker.AddTotalRecords(1)
	select {
	case s.inCh <- rec:
		return true
	case <-s.ctx.Done():
		s.tracker.AdjustTotalRecords(-1)
		s.budget.Release(1, size)
		return false
	}
}

func columnKind(typ columnv1.ColumnType) metrics.Kind {
	switch typ {
	case columnv1.ColumnType_COLUMN_TYPE_SCALAR:
		return metrics.KindScalar
	case columnv1.ColumnType_COLUMN_TYPE_HISTOGRAM:
		return metrics.KindHistogram
	default:
		return metrics.KindMedia
	}
}
