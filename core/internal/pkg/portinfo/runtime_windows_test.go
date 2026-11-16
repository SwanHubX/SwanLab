//go:build windows

package portinfo

import "testing"

func TestWindowsUserSegmentSanitizationAndEnvPriority(t *testing.T) {
	cases := []struct {
		name string
		env  map[string]string
		want string
	}{
		{
			name: "domain username",
			env:  map[string]string{"USERNAME": `CORP\Alice.Smith`},
			want: "CORP_Alice.Smith",
		},
		{
			name: "unicode only username",
			env:  map[string]string{"USERNAME": "张三"},
			want: "default",
		},
		{
			name: "mixed ascii and unicode",
			env:  map[string]string{"USERNAME": "alice_张三"},
			want: "alice",
		},
		{
			name: "dots and underscores only",
			env:  map[string]string{"USERNAME": "...___"},
			want: "default",
		},
		{
			name: "env fallback priority",
			env:  map[string]string{"USERNAME": "  ", "USER": "bob", "LOGNAME": "carol"},
			want: "bob",
		},
		{
			name: "all empty falls back to default",
			env:  map[string]string{},
			want: "default",
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := windowsUserSegment(func(k string) string { return tc.env[k] })
			if got != tc.want {
				t.Fatalf("windowsUserSegment = %q, want %q", got, tc.want)
			}
		})
	}
}
