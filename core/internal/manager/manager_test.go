package manager

import (
	"errors"
	"fmt"
	"sync"
	"testing"

	"github.com/google/uuid"
)

func TestManagerLifecycle(t *testing.T) {
	m := New()
	handle, err := m.Start()
	if err != nil || handle == "" {
		t.Fatalf("Start: handle=%q err=%v", handle, err)
	}
	if state, err := m.State(handle); err != nil || state != RunStateRunning {
		t.Fatalf("State after start: state=%v err=%v", state, err)
	}
	if err := m.Confirm(handle); !errors.Is(err, ErrRunNotFinished) {
		t.Fatalf("early Confirm err=%v, want ErrRunNotFinished", err)
	}
	if err := m.Finish(handle); err != nil {
		t.Fatalf("Finish: %v", err)
	}
	if state, err := m.State(handle); err != nil || state != RunStateFinished {
		t.Fatalf("State after finish: state=%v err=%v", state, err)
	}
	if err := m.Confirm(handle); err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if _, err := m.State(handle); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after confirm err=%v, want ErrSessionNotFound", err)
	}
}

func TestStartGeneratesUUIDv7(t *testing.T) {
	m := New()
	id, err := m.Start()
	if err != nil {
		t.Fatalf("Start: %v", err)
	}
	parsed, err := uuid.Parse(id)
	if err != nil {
		t.Fatalf("session id %q is not a uuid: %v", id, err)
	}
	if parsed.Version() != 7 {
		t.Fatalf("session id version = %d, want 7", parsed.Version())
	}
}

func TestFinishConfirmRace(t *testing.T) {
	m := New()
	for i := 0; i < 20; i++ {
		handle, err := m.Start()
		if err != nil {
			t.Fatal(err)
		}
		var wg sync.WaitGroup
		var finishErr, confirmErr error
		wg.Add(2)
		go func() {
			defer wg.Done()
			finishErr = m.Finish(handle)
		}()
		go func() {
			defer wg.Done()
			confirmErr = m.Confirm(handle)
		}()
		wg.Wait()
		if finishErr != nil {
			t.Fatalf("iteration %d Finish: %v", i, finishErr)
		}
		if confirmErr != nil && !errors.Is(confirmErr, ErrRunNotFinished) {
			t.Fatalf("iteration %d Confirm: %v", i, confirmErr)
		}
		if confirmErr != nil {
			if err := m.Confirm(handle); err != nil {
				t.Fatalf("iteration %d retry Confirm: %v", i, err)
			}
		}
	}
}

func TestRunIsolation(t *testing.T) {
	m := New()
	handleA, err := m.Start()
	if err != nil {
		t.Fatal(err)
	}
	handleB, err := m.Start()
	if err != nil {
		t.Fatal(err)
	}
	if handleA == handleB {
		t.Fatal("two runs must receive distinct handles")
	}
	if err := m.Finish(handleA); err != nil {
		t.Fatal(err)
	}
	if state, _ := m.State(handleA); state != RunStateFinished {
		t.Fatalf("run A state=%v, want Finished", state)
	}
	if state, _ := m.State(handleB); state != RunStateRunning {
		t.Fatalf("run B state=%v, want Running", state)
	}
}

func TestStartRejectsWhenFull(t *testing.T) {
	m := New()
	handles := make([]string, 0, maxSessions)
	for i := 0; i < maxSessions; i++ {
		handle, err := m.Start()
		if err != nil {
			t.Fatalf("start %d: %v", i, err)
		}
		handles = append(handles, handle)
	}
	if _, err := m.Start(); !errors.Is(err, ErrTooManySessions) {
		t.Fatalf("err = %v, want ErrTooManySessions", err)
	}
	// 释放一个会话后应重新可创建，且不影响其它会话
	if err := m.Finish(handles[0]); err != nil {
		t.Fatalf("Finish: %v", err)
	}
	if err := m.Confirm(handles[0]); err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if _, err := m.Start(); err != nil {
		t.Fatalf("Start after release: %v", err)
	}
	if state, err := m.State(handles[1]); err != nil || state != RunStateRunning {
		t.Fatalf("peer session state=%v err=%v, want Running", state, err)
	}
}

func TestConcurrentRuns(t *testing.T) {
	m := New()
	const clients = 4
	var wg sync.WaitGroup
	errs := make(chan error, clients)
	for i := 0; i < clients; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()
			handle, err := m.Start()
			if err != nil {
				errs <- fmt.Errorf("client %d start: %w", id, err)
				return
			}
			for j := 0; j < 20; j++ {
				if _, err = m.State(handle); err != nil {
					errs <- fmt.Errorf("client %d state: %w", id, err)
					return
				}
			}
			if err = m.Finish(handle); err != nil {
				errs <- fmt.Errorf("client %d finish: %w", id, err)
				return
			}
			if err = m.Confirm(handle); err != nil {
				errs <- fmt.Errorf("client %d confirm: %w", id, err)
			}
		}(i)
	}
	wg.Wait()
	close(errs)
	for err := range errs {
		t.Error(err)
	}
}
