package server

import (
	"crypto/rand"
	"encoding/base64"
	"sync"

	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/status"

	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
)

// runSessionState 是单个 run 会话的内部状态（不上 proto、不并入 CoreState）：
// running → finishStaged（DeliverRunFinish）→ confirmed/释放（ConfirmRunFinish）。
type runSessionState int32

const (
	sessionRunning runSessionState = iota
	sessionFinishStaged
)

// runSession 记录单个 run 会话的路由状态。
type runSession struct {
	mu    sync.Mutex
	state runSessionState
}

// finish 幂等标记 finish 已交付。
func (s *runSession) finish() {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.state = sessionFinishStaged
}

// finished 报告 finish 是否完成 record 上传
// ConfirmRunFinish 的前置条件
func (s *runSession) finished() bool {
	s.mu.Lock()
	defer s.mu.Unlock()
	return s.state == sessionFinishStaged
}

// coreState 把内部会话状态映射为线上契约 CoreState（run 级数据排空轴）。
// finish 交付后骨架无数据可排空，直接报 FINISHED，drain 轮询即可判停。
func (s *runSession) coreState() operationv1.CoreState {
	s.mu.Lock()
	defer s.mu.Unlock()
	if s.state == sessionRunning {
		return operationv1.CoreState_CORE_STATE_RUNNING
	}
	return operationv1.CoreState_CORE_STATE_FINISHED
}

// sessionRegistry 是并发安全的 run_handle → 会话映射。
// DeliverRunStart 创建，ConfirmRunFinish 释放；confirm 绝不关闭 gRPC Server。
type sessionRegistry struct {
	mu       sync.RWMutex
	sessions map[string]*runSession
}

func newSessionRegistry() *sessionRegistry {
	return &sessionRegistry{sessions: make(map[string]*runSession)}
}

// lookup 校验并返回目标会话；空 handle 为 InvalidArgument，
// 未知或已释放的 handle 为 NotFound。
func (r *sessionRegistry) lookup(handle string) (*runSession, error) {
	if handle == "" {
		return nil, status.Error(codes.InvalidArgument, "run_handle must not be empty")
	}
	r.mu.RLock()
	defer r.mu.RUnlock()
	s, ok := r.sessions[handle]
	if !ok {
		return nil, status.Error(codes.NotFound, "unknown or released run_handle")
	}
	return s, nil
}

// create 生成随机 opaque handle 并登记会话。
func (r *sessionRegistry) create() (string, *runSession, error) {
	buf := make([]byte, 32)
	if _, err := rand.Read(buf); err != nil {
		return "", nil, status.Error(codes.Internal, "generate run handle")
	}
	handle := base64.RawURLEncoding.EncodeToString(buf)
	s := &runSession{}
	r.mu.Lock()
	defer r.mu.Unlock()
	if _, exists := r.sessions[handle]; exists {
		// 256-bit 随机碰撞概率可忽略；出现即视为内部错误。
		return "", nil, status.Error(codes.Internal, "run handle collision")
	}
	r.sessions[handle] = s
	return handle, s, nil
}

// confirm 在 registry 锁下原子校验并摘除会话：handle 为空返回
// InvalidArgument，未知或已释放返回 NotFound，finish 未交付返回
// FailedPrecondition 且会话保留。confirm 成功时 finish 已生效。
// 锁顺序：registry mu → session mu。
func (r *sessionRegistry) confirm(handle string) error {
	if handle == "" {
		return status.Error(codes.InvalidArgument, "run_handle must not be empty")
	}
	r.mu.Lock()
	defer r.mu.Unlock()
	s, ok := r.sessions[handle]
	if !ok {
		return status.Error(codes.NotFound, "unknown or released run_handle")
	}
	if !s.finished() {
		return status.Error(codes.FailedPrecondition, "run finish has not been delivered for this handle")
	}
	delete(r.sessions, handle)
	return nil
}
