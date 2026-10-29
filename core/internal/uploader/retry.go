package uploader

import (
	"math"
	"math/rand"
	"time"
)

// 默认退避参数，对齐 Python RetryPolicy。
const (
	// DefaultRetryBase 是首次失败后的等待时长。
	DefaultRetryBase = time.Second
	// DefaultRetryMultiplier 是每轮失败的放大倍数。
	DefaultRetryMultiplier = 2.0
	// DefaultRetryMax 是单次等待的上限。
	DefaultRetryMax = 10 * time.Second
	// DefaultRetryJitterLow 是抖动下界。
	DefaultRetryJitterLow = 0.8
	// DefaultRetryJitterHigh 是抖动上界。
	DefaultRetryJitterHigh = 1.2
)

// RetryPolicy 计算单轨连续失败后的等待时长（指数退避 + 抖动），可并发使用。
type RetryPolicy struct {
	Base       time.Duration  // 首轮等待，<=0 取 DefaultRetryBase
	Multiplier float64        // 每轮放大倍数，<=0 取 DefaultRetryMultiplier
	Max        time.Duration  // 单次等待上限，<=0 取 DefaultRetryMax
	JitterLow  float64        // 抖动下界，JitterHigh<=0 时连同 JitterHigh 一起取默认值
	JitterHigh float64        // 抖动上界
	Rand       func() float64 // 返回 [0,1] 的随机源，nil 取 rand.Float64
}

// DefaultRetryPolicy 返回与 Python 一致的退避参数。
func DefaultRetryPolicy() RetryPolicy {
	return RetryPolicy{
		Base:       DefaultRetryBase,
		Multiplier: DefaultRetryMultiplier,
		Max:        DefaultRetryMax,
		JitterLow:  DefaultRetryJitterLow,
		JitterHigh: DefaultRetryJitterHigh,
		Rand:       rand.Float64,
	}
}

// Delay 返回第 failedRounds 轮失败后的等待时长（<1 按 1 处理）。
// 指数项先饱和并截断到 Max，再乘抖动后二次截断，轮次极大也不会溢出。
func (p RetryPolicy) Delay(failedRounds int) time.Duration {
	p = p.withDefaults()
	exponent := min(max(failedRounds, 1)-1, p.saturationExponent())

	delay := float64(p.Base) * math.Pow(p.Multiplier, float64(exponent))
	delay = math.Min(delay, float64(p.Max))
	delay *= p.jitter()
	return time.Duration(math.Round(math.Min(delay, float64(p.Max))))
}

// isZero 报告是否为零值。
func (p RetryPolicy) isZero() bool {
	return p.Base == 0 && p.Multiplier == 0 && p.Max == 0 &&
		p.JitterLow == 0 && p.JitterHigh == 0 && p.Rand == nil
}

// withDefaults 把零值字段补成默认值。
func (p RetryPolicy) withDefaults() RetryPolicy {
	if p.Base <= 0 {
		p.Base = DefaultRetryBase
	}
	if p.Multiplier <= 0 {
		p.Multiplier = DefaultRetryMultiplier
	}
	if p.Max <= 0 {
		p.Max = DefaultRetryMax
	}
	if p.JitterHigh <= 0 {
		p.JitterLow = DefaultRetryJitterLow
		p.JitterHigh = DefaultRetryJitterHigh
	} else if p.JitterLow < 0 {
		p.JitterLow = 0
	}
	if p.Rand == nil {
		p.Rand = rand.Float64
	}
	return p
}

// saturationExponent 是达到 Max 所需的指数上限。
func (p RetryPolicy) saturationExponent() int {
	if p.Multiplier <= 1 || p.Max <= p.Base {
		return 0
	}
	return int(math.Ceil(math.Log(float64(p.Max)/float64(p.Base)) / math.Log(p.Multiplier)))
}

// jitter 返回抖动区间内的均匀系数。
func (p RetryPolicy) jitter() float64 {
	r := math.Min(math.Max(p.Rand(), 0), 1)
	return p.JitterLow + r*(p.JitterHigh-p.JitterLow)
}
