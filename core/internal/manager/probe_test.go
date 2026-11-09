package manager

import (
	"context"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
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
		if got := probeReachable(context.Background(), srv.URL+c.path); got != c.want {
			t.Fatalf("probeReachable(%s) = %v, want %v", c.path, got, c.want)
		}
	}

	closed := httptest.NewServer(http.HandlerFunc(func(http.ResponseWriter, *http.Request) {}))
	url := closed.URL
	closed.Close()
	if probeReachable(context.Background(), url) {
		t.Fatal("closed server must be unreachable")
	}
}

func TestSiteKeyHidesAPIKey(t *testing.T) {
	a := siteKey("https://api.example.com", "secret-key")
	b := siteKey("https://api.example.com", "other-key")
	if a == b {
		t.Fatal("different api keys must not share a site")
	}
	if len(a) != 32 || strings.Contains(a, "secret-key") {
		t.Fatalf("site key %q leaks the api key or has unexpected length", a)
	}
}
