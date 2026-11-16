//go:build windows

package portinfo

import (
	"os"
	"os/user"
	"strings"
)

func defaultRuntimeRoot() string {
	for _, key := range []string{"TEMP", "TMP"} {
		if value := strings.TrimSpace(os.Getenv(key)); value != "" {
			return value
		}
	}
	return os.TempDir()
}

func currentUserSegment() (string, error) {
	if current, err := user.Current(); err == nil && strings.TrimSpace(current.Username) != "" {
		return sanitizeUserSegment(current.Username), nil
	}
	if username := strings.TrimSpace(os.Getenv("USERNAME")); username != "" {
		return sanitizeUserSegment(username), nil
	}
	return "default", nil
}

func sanitizeUserSegment(raw string) string {
	var b strings.Builder
	for _, r := range raw {
		switch {
		case (r >= 'a' && r <= 'z') || (r >= 'A' && r <= 'Z') || (r >= '0' && r <= '9') || r == '-' || r == '_' || r == '.':
			b.WriteRune(r)
		default:
			b.WriteByte('_')
		}
	}
	cleaned := strings.Trim(b.String(), "._")
	if cleaned == "" {
		return "default"
	}
	return cleaned
}

func ensureRuntimeDirPerm(string) error {
	return nil
}

func validateRuntimeDirOwner(string, os.FileInfo) error {
	return nil
}

func validateSocketPathPlatform(string) error {
	return nil
}
