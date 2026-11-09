package health

import (
	"sync"
	"time"
)

// SiteOptions 描述一个站点的探针配置；同一 Key 的后续 Acquire 复用首个实例。
type SiteOptions struct {
	// Key 是站点标识，约定为 api_host + api_key 摘要。
	Key string
	// Interval / Timeout 是探测周期与单次探测超时，<=0 取默认值。
	Interval time.Duration
	Timeout  time.Duration
	// DegradeThreshold / RecoverThreshold 是连续失败 / 成功阈值，<=0 取默认值。
	DegradeThreshold int
	RecoverThreshold int
	// Control / Data 分别是控制面与数据面的探测函数，必须非 nil。
	Control ProbeFunc
	Data    ProbeFunc
}

// Registry 按站点共享探针：引用计数归零时停止探测协程。
type Registry struct {
	clock Clock

	mu    sync.Mutex
	sites map[string]*entry
}

type entry struct {
	site *Site
	refs int
}

// NewRegistry 创建注册表；clock 为 nil 时使用真实定时器。
func NewRegistry(clock Clock) *Registry {
	if clock == nil {
		clock = realClock{}
	}
	return &Registry{clock: clock, sites: make(map[string]*entry)}
}

// Acquire 增加站点的引用计数并返回共享实例；首次获取时启动两个探测协程。
// Key 为空或任一探针为 nil 属于调用方错误，直接 panic。
func (r *Registry) Acquire(opts SiteOptions) *Site {
	if opts.Key == "" || opts.Control == nil || opts.Data == nil {
		panic("health: Acquire requires key and both probe functions")
	}
	r.mu.Lock()
	defer r.mu.Unlock()
	if e, ok := r.sites[opts.Key]; ok {
		e.refs++
		return e.site
	}
	s := newSite(opts, r.clock)
	r.sites[opts.Key] = &entry{site: s, refs: 1}
	return s
}

// Release 减少引用计数；归零时停止探针并移除站点。nil 或重复归还可安全调用。
func (r *Registry) Release(s *Site) {
	if s == nil {
		return
	}
	r.mu.Lock()
	defer r.mu.Unlock()
	e, ok := r.sites[s.key]
	if !ok || e.site != s {
		return
	}
	e.refs--
	if e.refs > 0 {
		return
	}
	delete(r.sites, s.key)
	e.site.stop()
}

// Len 返回当前存活的站点数，共享站点只计一次。
func (r *Registry) Len() int {
	r.mu.Lock()
	defer r.mu.Unlock()
	return len(r.sites)
}
