package config

import (
	"errors"
	"math"
	"strings"
	"testing"
	"time"
)

func TestDefaults(t *testing.T) {
	s := Defaults()
	if s.RecordBatch != DefaultRecordBatch {
		t.Fatalf("RecordBatch = %d, want %d", s.RecordBatch, DefaultRecordBatch)
	}
	if s.RecordInterval != DefaultRecordInterval {
		t.Fatalf("RecordInterval = %v, want %v", s.RecordInterval, DefaultRecordInterval)
	}
	if s.SaveSize != DefaultSaveSize || s.SaveSplit != DefaultSaveSplit ||
		s.SavePart != DefaultSavePart || s.SaveBatch != DefaultSaveBatch {
		t.Fatalf("save defaults mismatch: %+v", s)
	}
	if s.UploadBudgetRecords != DefaultUploadBudgetRecords || s.UploadBudgetBytes != DefaultUploadBudgetBytes {
		t.Fatalf("upload budget defaults = (%d, %d), want (%d, %d)",
			s.UploadBudgetRecords, s.UploadBudgetBytes, DefaultUploadBudgetRecords, DefaultUploadBudgetBytes)
	}
	// 这些字段默认为零值。
	if s.RunID != "" || s.RunDir != "" || s.SectionRule != 0 || s.SkipStore || s.RecoverOnline || s.APIKey != "" {
		t.Fatalf("unexpected non-zero defaults: %+v", s)
	}
	// Defaults 必须能通过校验（run_id / run_dir 由调用方补齐）。
	s.RunID, s.RunDir = "r", "/rd"
	if err := s.Validate(); err != nil {
		t.Fatalf("Defaults + 必填项应通过校验: %v", err)
	}
}

func TestValidate(t *testing.T) {
	base := func() Settings {
		s := Defaults()
		s.RunID, s.RunDir = "r", "/rd"
		return s
	}

	tests := []struct {
		name   string
		mutate func(*Settings)
		wantIs error
		want   string
	}{
		{"ok", func(*Settings) {}, nil, ""},
		{"缺少 run_id", func(s *Settings) { s.RunID = "" }, ErrRunIDRequired, ""},
		{"缺少 run_dir", func(s *Settings) { s.RunDir = "" }, ErrRunDirRequired, ""},
		{"record_batch 为 0", func(s *Settings) { s.RecordBatch = 0 }, nil, "record_batch"},
		{"record_batch 为负", func(s *Settings) { s.RecordBatch = -1 }, nil, "record_batch"},
		{"record_batch 达上界", func(s *Settings) { s.RecordBatch = recordBatchMax }, nil, "record_batch"},
		{"record_batch 上界内", func(s *Settings) { s.RecordBatch = recordBatchMax - 1 }, nil, ""},
		{"record_interval 为 0", func(s *Settings) { s.RecordInterval = 0 }, nil, "record_interval"},
		{"record_interval 为负", func(s *Settings) { s.RecordInterval = -1.5 }, nil, "record_interval"},
		{"record_interval 为 NaN", func(s *Settings) { s.RecordInterval = math.NaN() }, nil, "record_interval"},
		{"record_interval 为正无穷", func(s *Settings) { s.RecordInterval = math.Inf(1) }, nil, ""},
		{"record_interval 为负无穷", func(s *Settings) { s.RecordInterval = math.Inf(-1) }, nil, "record_interval"},
		{"save_part 为 0", func(s *Settings) { s.SavePart = 0 }, nil, "save_part"},
		// save_size / save_split / save_batch 不额外约束。
		{"save_batch 无 bounds", func(s *Settings) { s.SaveBatch = -1 }, nil, ""},
		{"save_size 无 bounds", func(s *Settings) { s.SaveSize = -1 }, nil, ""},
		{"save_split 无 bounds", func(s *Settings) { s.SaveSplit = -1 }, nil, ""},
		{"section_rule 允许负值", func(s *Settings) { s.SectionRule = -3 }, nil, ""},
		{"mode 未知", func(s *Settings) { s.Mode = "weird" }, nil, "mode"},
		{"local 允许", func(s *Settings) { s.Mode = ModeLocal }, nil, ""},
		{"local 禁止 skip_store", func(s *Settings) { s.Mode = ModeLocal; s.SkipStore = true }, nil, "skip_store"},
		{"offline 禁止 skip_store", func(s *Settings) { s.Mode = ModeOffline; s.SkipStore = true }, nil, "skip_store"},
		{"disabled 允许", func(s *Settings) { s.Mode = ModeDisabled }, nil, ""},
		{"online 允许 recover_online", func(s *Settings) { s.Mode = ModeOnline; s.RecoverOnline = true }, nil, ""},
		{"空 mode 允许 recover_online", func(s *Settings) { s.Mode = ""; s.RecoverOnline = true }, nil, ""},
		{"skip_store 禁止 recover_online", func(s *Settings) { s.Mode = ModeOnline; s.SkipStore = true; s.RecoverOnline = true }, nil, "recover_online"},
		{"local 禁止 recover_online", func(s *Settings) { s.Mode = ModeLocal; s.RecoverOnline = true }, nil, "recover_online"},
		{"offline 禁止 recover_online", func(s *Settings) { s.Mode = ModeOffline; s.RecoverOnline = true }, nil, "recover_online"},
		{"disabled 禁止 recover_online", func(s *Settings) { s.Mode = ModeDisabled; s.RecoverOnline = true }, nil, "recover_online"},
		{"health_probe_interval 为 0", func(s *Settings) { s.HealthProbeInterval = 0 }, nil, ""},
		{"health_probe_interval 为负", func(s *Settings) { s.HealthProbeInterval = -0.1 }, nil, "health_probe_interval"},
		{"health_probe_interval 为 NaN", func(s *Settings) { s.HealthProbeInterval = math.NaN() }, nil, "health_probe_interval"},
		{"health_probe_interval 为正无穷", func(s *Settings) { s.HealthProbeInterval = math.Inf(1) }, nil, "health_probe_interval"},
		{"health_probe_timeout 为负", func(s *Settings) { s.HealthProbeTimeout = -1 }, nil, "health_probe_timeout"},
		{"health_probe_timeout 为 NaN", func(s *Settings) { s.HealthProbeTimeout = math.NaN() }, nil, "health_probe_timeout"},
		{"health_probe_timeout 为正无穷", func(s *Settings) { s.HealthProbeTimeout = math.Inf(1) }, nil, "health_probe_timeout"},
		{"health_degrade_threshold 为负", func(s *Settings) { s.HealthDegradeThreshold = -1 }, nil, "health thresholds"},
		{"health_recover_threshold 为负", func(s *Settings) { s.HealthRecoverThreshold = -1 }, nil, "health thresholds"},
		{"health 参数合法", func(s *Settings) {
			s.HealthProbeInterval = 1.5
			s.HealthProbeTimeout = 0.75
			s.HealthDegradeThreshold = 5
			s.HealthRecoverThreshold = 4
		}, nil, ""},
		{"budget 允许非正", func(s *Settings) {
			s.UploadBudgetRecords = -1
			s.UploadBudgetBytes = -1
		}, nil, ""},
		{"budget 允许正记录数", func(s *Settings) { s.UploadBudgetRecords = 1 }, nil, ""},
		{"budget 允许正字节数", func(s *Settings) { s.UploadBudgetBytes = 1 }, nil, ""},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			s := base()
			tc.mutate(&s)
			err := s.Validate()
			if tc.wantIs != nil {
				if !errors.Is(err, tc.wantIs) {
					t.Fatalf("Validate() = %v, want errors.Is %v", err, tc.wantIs)
				}
				return
			}
			if tc.want == "" {
				if err != nil {
					t.Fatalf("Validate() = %v, want nil", err)
				}
				return
			}
			if err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Fatalf("Validate() = %v, want 包含 %q", err, tc.want)
			}
		})
	}
}

func TestNormalizedMode(t *testing.T) {
	tests := []struct {
		mode Mode
		want Mode
	}{
		{"", ModeOnline},
		{ModeOnline, ModeOnline},
		{ModeLocal, ModeLocal},
		{ModeOffline, ModeOffline},
		{ModeDisabled, ModeDisabled},
	}
	for _, tc := range tests {
		if got := (Settings{Mode: tc.mode}).NormalizedMode(); got != tc.want {
			t.Fatalf("NormalizedMode(%q) = %q, want %q", tc.mode, got, tc.want)
		}
	}
}

func TestRecordIntervalDuration(t *testing.T) {
	tests := []struct {
		name     string
		interval float64
		want     time.Duration
	}{
		{"默认间隔", DefaultRecordInterval, 1500 * time.Millisecond},
		{"亚毫秒抬到 1ms", 0.0001, time.Millisecond},
		{"普通值", 0.25, 250 * time.Millisecond},
		// +Inf 与超大值不能溢出成负数（否则 uploader 会当成未设置而静默回退 1.5s）。
		{"正无穷截断到上限", math.Inf(1), maxRecordInterval},
		{"超大值截断到上限", 1e18, maxRecordInterval},
		{"恰好上限", maxRecordInterval.Seconds(), maxRecordInterval},
		{"零值取默认", 0, 1500 * time.Millisecond},
		{"负值取默认", -1, 1500 * time.Millisecond},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := (Settings{RecordInterval: tc.interval}).RecordIntervalDuration()
			if got != tc.want {
				t.Fatalf("RecordIntervalDuration(%v) = %v, want %v", tc.interval, got, tc.want)
			}
			if got <= 0 {
				t.Fatalf("RecordIntervalDuration(%v) = %v, 不应当为负数", tc.interval, got)
			}
		})
	}
}

func TestMaskKey(t *testing.T) {
	tests := []struct {
		in, want string
	}{
		{"", ""},
		{"short", "****"},
		{"12345678", "****"},
		{"123456789", "1234...6789"},
		{"swanlab-api-key-0123456789", "swan...6789"},
	}
	for _, tc := range tests {
		if got := MaskKey(tc.in); got != tc.want {
			t.Fatalf("MaskKey(%q) = %q, want %q", tc.in, got, tc.want)
		}
	}
}

func TestStringRedactsAPIKey(t *testing.T) {
	s := Defaults()
	s.RunID, s.RunDir = "r", "/rd"
	s.APIKey = "super-secret-key-000000"
	s.APIHost = "https://api.example.com"
	s.WebHost = "https://web.example.com"

	out := s.String()
	if strings.Contains(out, s.APIKey) {
		t.Fatalf("String() 泄漏了完整 API Key: %s", out)
	}
	if !strings.Contains(out, MaskKey(s.APIKey)) {
		t.Fatalf("String() 未包含脱敏 key: %s", out)
	}
	for _, want := range []string{s.RunID, s.APIHost, s.WebHost, "record_batch", "save_part", "recover_online=false"} {
		if !strings.Contains(out, want) {
			t.Fatalf("String() 缺少字段 %q: %s", want, out)
		}
	}
	s.RecoverOnline = true
	if outRecover := s.String(); !strings.Contains(outRecover, "recover_online=true") {
		t.Fatalf("String() 未包含 recover_online=true: %s", outRecover)
	}
}
