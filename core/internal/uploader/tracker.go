package uploader

import (
	"math"
	"sort"
	"sync"
	"time"

	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
)

const (
	// trackerRateAlpha 是 EMA 平滑因子，越大越贴近近期吞吐。
	trackerRateAlpha = 0.3
	// trackerRateWindow 是速率采样窗口：窗口内的增量累积后一次性混合。
	trackerRateWindow = 500 * time.Millisecond
	// defaultVisibleFiles 是快照里展示的文件条数上限。
	defaultVisibleFiles = 5
)

// rateTracker 是 EMA 速率追踪器，对齐 Python transport/tracker/_RateTracker。
//
// 窗口内的增量先累积，到窗口边界才换算成瞬时速率并混合，避免逐字节抖动；
// 快照读取时按距上次混合的间隔做指数衰减，长时间无进展的速率会自然归零。
type rateTracker struct {
	rate    float64
	last    time.Time
	pending int64
}

// advance 累积 delta，并在距上次混合满一个窗口时更新速率。
func (r *rateTracker) advance(now time.Time, delta int64) {
	if delta < 0 {
		return
	}
	r.pending += delta
	if r.last.IsZero() {
		r.last = now
		return
	}
	elapsed := now.Sub(r.last)
	if elapsed < trackerRateWindow {
		return
	}
	instant := float64(r.pending) / elapsed.Seconds()
	if r.rate <= 0 {
		r.rate = instant
	} else {
		r.rate = trackerRateAlpha*instant + (1-trackerRateAlpha)*r.rate
	}
	r.last = now
	r.pending = 0
}

// decayed 返回经指数衰减后的当前速率。
func (r *rateTracker) decayed(now time.Time) float64 {
	if r.rate <= 0 || r.last.IsZero() {
		return 0
	}
	elapsed := now.Sub(r.last)
	if elapsed < trackerRateWindow {
		return r.rate
	}
	intervals := elapsed.Seconds() / trackerRateWindow.Seconds()
	return r.rate * math.Pow(1-trackerRateAlpha, intervals)
}

// touch 在首次活动时记录起点，使首个窗口有可用的时间基准。
func (r *rateTracker) touch(now time.Time) {
	if r.last.IsZero() {
		r.last = now
	}
}

// trackedFile 是单个文件的进度，uploaded 由分片增量累计得到。
type trackedFile struct {
	path     string
	total    int64
	uploaded int64
	sequence int
	parts    map[string]int64
	rate     rateTracker
}

// Tracker 是并发安全的上传进度账本，对齐 Python transport/tracker/UploadTracker。
//
// 记录维度供上传管线记账（总记录数 / 已终结记录数），文件维度供 PR-6 使用；
// Snapshot 直接产出 operationv1.OperationStats，供心跳与状态查询共用。
// 零值不可用，应由 NewTracker 构造。
type Tracker struct {
	mu sync.Mutex

	now   func() time.Time
	limit int

	state           operationv1.CoreState
	totalNumber     int64
	uploadedNumber  int64
	totalSize       int64
	uploadedSize    int64
	totalRecords    int64
	uploadedRecords int64

	globalRate rateTracker
	files      map[string]*trackedFile
	registered map[string]struct{}
	sequence   int
}

// NewTracker 构造默认追踪器：可见文件上限 5，时钟取 time.Now。
func NewTracker() *Tracker {
	return &Tracker{
		now:        time.Now,
		limit:      defaultVisibleFiles,
		files:      make(map[string]*trackedFile),
		registered: make(map[string]struct{}),
	}
}

// SetState 设置 Core 生命周期状态。
func (t *Tracker) SetState(state operationv1.CoreState) {
	t.mu.Lock()
	defer t.mu.Unlock()
	t.state = state
}

// AddTotalRecords 累加已知总记录数；n<=0 忽略。
func (t *Tracker) AddTotalRecords(n int64) {
	if n <= 0 {
		return
	}
	t.mu.Lock()
	defer t.mu.Unlock()
	t.totalRecords += n
}

// AdvanceRecords 累加已终结记录数；n<=0 忽略。
func (t *Tracker) AdvanceRecords(n int64) {
	if n <= 0 {
		return
	}
	t.mu.Lock()
	defer t.mu.Unlock()
	t.uploadedRecords += n
}

// AddFile 注册一个待上传文件；key 幂等，total<=0 忽略。
func (t *Tracker) AddFile(key, path string, total int64) {
	if total <= 0 {
		return
	}
	t.mu.Lock()
	defer t.mu.Unlock()
	if _, ok := t.registered[key]; ok {
		return
	}
	now := t.now()
	t.registered[key] = struct{}{}
	t.sequence++
	t.files[key] = &trackedFile{
		path:     path,
		total:    total,
		sequence: t.sequence,
		parts:    make(map[string]int64),
		rate:     rateTracker{last: now},
	}
	t.totalNumber++
	t.totalSize += total
	t.globalRate.touch(now)
}

// UpdateFileProgress 上报某个分片的上传进度，按分片 key 计算增量。
// 负增量且 current>0 视为乱序上报并忽略；current==0 允许把分片重置。
func (t *Tracker) UpdateFileProgress(fileKey, uploadKey string, current int64) {
	current = normalizeAmount(current)
	t.mu.Lock()
	defer t.mu.Unlock()
	file, ok := t.files[fileKey]
	if !ok {
		return
	}
	prev := file.parts[uploadKey]
	delta := current - prev
	if delta == 0 || (delta < 0 && current > 0) {
		return
	}
	file.parts[uploadKey] = current
	file.uploaded += delta

	now := t.now()
	file.rate.advance(now, delta)
	t.uploadedSize += delta
	t.globalRate.advance(now, delta)
}

// FinishFile 标记文件上传完成：补齐未上报的剩余字节，并计入完成文件数。
func (t *Tracker) FinishFile(key string) {
	t.mu.Lock()
	defer t.mu.Unlock()
	file, ok := t.files[key]
	if !ok {
		return
	}
	delete(t.files, key)
	if remaining := file.total - file.uploaded; remaining > 0 {
		t.uploadedSize += remaining
	}
	t.uploadedNumber++
}

// FailFile 移除已终结但未成功上传的文件，不计入完成文件数与字节。
func (t *Tracker) FailFile(key string) {
	t.mu.Lock()
	defer t.mu.Unlock()
	delete(t.files, key)
}

// Snapshot 返回当前进度快照；返回的 OperationStats 归调用方所有。
func (t *Tracker) Snapshot() *operationv1.OperationStats {
	t.mu.Lock()
	defer t.mu.Unlock()

	now := t.now()
	stats := &operationv1.OperationStats{
		State:           t.state,
		TotalNumber:     t.totalNumber,
		UploadedNumber:  t.uploadedNumber,
		TotalRecords:    t.totalRecords,
		UploadedRecords: t.uploadedRecords,
		TotalSize:       t.totalSize,
		UploadedSize:    t.uploadedSize,
		Rate:            t.globalRate.decayed(now),
	}
	visible := t.visibleFiles()
	if len(visible) == 0 {
		return stats
	}
	stats.Files = make([]*operationv1.FileProgress, 0, len(visible))
	for _, file := range visible {
		stats.Files = append(stats.Files, &operationv1.FileProgress{
			Path:     file.path,
			Total:    file.total,
			Uploaded: file.uploaded,
			Rate:     file.rate.decayed(now),
		})
	}
	return stats
}

// visibleFiles 返回未完成文件中前 limit 个：已上传过字节的排在前面，同组按注册序。
func (t *Tracker) visibleFiles() []*trackedFile {
	if t.limit <= 0 {
		return nil
	}
	visible := make([]*trackedFile, 0, len(t.files))
	for _, file := range t.files {
		if file.uploaded < file.total {
			visible = append(visible, file)
		}
	}
	sort.Slice(visible, func(i, j int) bool {
		if (visible[i].uploaded <= 0) != (visible[j].uploaded <= 0) {
			return visible[i].uploaded > 0
		}
		return visible[i].sequence < visible[j].sequence
	})
	if len(visible) > t.limit {
		visible = visible[:t.limit]
	}
	return visible
}

// normalizeAmount 把非正数归零，保证计数只增不减。
func normalizeAmount(n int64) int64 {
	if n <= 0 {
		return 0
	}
	return n
}
