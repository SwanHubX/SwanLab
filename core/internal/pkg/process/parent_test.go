package process

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

func TestNotifyOnParentExitRejectsInvalidPID(t *testing.T) {
	console.Reset()
	defer console.Reset()
	dir := t.TempDir()
	if err := console.Init(dir); err != nil {
		t.Fatalf("console.Init: %v", err)
	}

	for _, parentPID := range []int{-1, 0} {
		parentExited, err := NotifyOnParentExit(parentPID)
		if err == nil {
			t.Fatalf("NotifyOnParentExit(%d) returned nil error", parentPID)
		}
		if parentExited != nil {
			t.Fatalf("NotifyOnParentExit(%d) channel is non-nil on error", parentPID)
		}
	}

	data, err := os.ReadFile(filepath.Join(dir, "debug-core.log"))
	if err != nil {
		t.Fatalf("read console log: %v", err)
	}
	if count := strings.Count(string(data), "Failed to monitor parent process: parent PID must be positive"); count != 2 {
		t.Fatalf("console log contains %d setup errors, want 2: %q", count, data)
	}
}

func TestNotifyOnParentExitDetectsChangedParent(t *testing.T) {
	parentExited, err := NotifyOnParentExit(os.Getppid() + 1)
	if !errors.Is(err, ErrParentExited) {
		t.Fatalf("NotifyOnParentExit() error = %v, want ErrParentExited", err)
	}
	if parentExited != nil {
		t.Fatal("NotifyOnParentExit() channel is non-nil on error")
	}
}
