package uploader

import (
	"testing"

	"go.uber.org/goleak"
)

// TestMain 断言 uploader 自己的协程全部退出；net/http 连接池的常驻协程不计入。
func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m,
		goleak.IgnoreTopFunction("net/http.(*persistConn).readLoop"),
		goleak.IgnoreTopFunction("net/http.(*persistConn).writeLoop"),
	)
}
