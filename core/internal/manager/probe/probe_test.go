package probe

import (
	"context"
	"io"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
	"time"
)

func TestProbeReachableClassification(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		switch r.URL.Path {
		case "/ok":
			w.WriteHeader(http.StatusNoContent)
		case "/missing":
			w.WriteHeader(http.StatusNotFound)
		case "/server-error":
			w.WriteHeader(http.StatusServiceUnavailable)
		case "/rate-limited":
			w.WriteHeader(http.StatusTooManyRequests)
		default:
			w.WriteHeader(http.StatusOK)
		}
	}))
	defer srv.Close()

	cases := []struct {
		path string
		want bool
	}{
		{"/ok", true},
		{"/missing", true},
		{"/server-error", false},
		{"/rate-limited", false},
	}
	for _, c := range cases {
		if got := Reachable(context.Background(), srv.URL+c.path); got != c.want {
			t.Fatalf("Reachable(%s) = %v, want %v", c.path, got, c.want)
		}
	}

	closed := httptest.NewServer(http.HandlerFunc(func(http.ResponseWriter, *http.Request) {}))
	targetURL := closed.URL
	closed.Close()
	if Reachable(context.Background(), targetURL) {
		t.Fatal("closed server must be unreachable")
	}
}

type roundTripFunc func(*http.Request) (*http.Response, error)

func (f roundTripFunc) RoundTrip(r *http.Request) (*http.Response, error) { return f(r) }

func TestProbeReachableHonorsConfiguredTimeout(t *testing.T) {
	old := probeClient.Transport
	t.Cleanup(func() { probeClient.Transport = old })

	var remaining time.Duration
	probeClient.Transport = roundTripFunc(func(r *http.Request) (*http.Response, error) {
		deadline, ok := r.Context().Deadline()
		if !ok {
			t.Fatal("probe request has no context deadline")
		}
		remaining = time.Until(deadline)
		return &http.Response{
			StatusCode: http.StatusOK,
			Body:       io.NopCloser(strings.NewReader("")),
			Header:     make(http.Header),
		}, nil
	})

	// 10s 预算必须原样传到请求；旧实现会被 Client.Timeout=2s 截断。
	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	defer cancel()
	if !Reachable(ctx, "http://probe.invalid") {
		t.Fatal("probe should be reachable")
	}
	if remaining < 5*time.Second {
		t.Fatalf("request deadline %v, want > 5s", remaining)
	}
}

func TestSiteKeyHidesAPIKey(t *testing.T) {
	a := SiteKey("https://api.example.com", "secret-key")
	b := SiteKey("https://api.example.com", "other-key")
	if a == b {
		t.Fatal("different api keys must not share a site")
	}
	if len(a) != 32 || strings.Contains(a, "secret-key") {
		t.Fatalf("site key %q leaks the api key or has unexpected length", a)
	}
}

func TestSecondsDuration(t *testing.T) {
	cases := []struct {
		in   float64
		want time.Duration
	}{
		{0, 0},
		{-1, 0},
		{1.5, 1500 * time.Millisecond},
		{0.75, 750 * time.Millisecond},
	}
	for _, c := range cases {
		if got := SecondsDuration(c.in); got != c.want {
			t.Fatalf("SecondsDuration(%v) = %v, want %v", c.in, got, c.want)
		}
	}
}
