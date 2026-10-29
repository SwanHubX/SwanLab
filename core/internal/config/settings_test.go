package config

import (
	"errors"
	"math"
	"strings"
	"testing"
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
	// 这些字段默认为零值。
	if s.RunID != "" || s.RunDir != "" || s.SectionRule != 0 || s.SkipStore || s.APIKey != "" {
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
	for _, want := range []string{s.RunID, s.APIHost, s.WebHost, "record_batch", "save_part"} {
		if !strings.Contains(out, want) {
			t.Fatalf("String() 缺少字段 %q: %s", want, out)
		}
	}
}
