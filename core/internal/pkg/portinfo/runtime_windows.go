//go:build windows

package portinfo

import (
	"os"
	"strings"
)

var windowsUserEnvKeys = []string{"USERNAME", "USER", "LOGNAME", "LNAME"}

func defaultRuntimeRoot() string {
	for _, key := range []string{"TEMP", "TMP"} {
		if value := strings.TrimSpace(os.Getenv(key)); value != "" {
			return value
		}
	}
	return os.TempDir()
}

func currentUserSegment() (string, error) {
	return windowsUserSegment(os.Getenv), nil
}

func windowsUserSegment(getenv func(string) string) string {
	for _, key := range windowsUserEnvKeys {
		if username := strings.TrimSpace(getenv(key)); username != "" {
			return sanitizeUserSegment(username)
		}
	}
	return "default"
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
