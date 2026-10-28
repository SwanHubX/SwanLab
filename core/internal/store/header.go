// Package store 实现 SwanLab 运行日志（run-*.swanlab）的本地持久化读写。
//
// 文件格式与 Python 侧 swanlab/sdk/internal/core_python/store 完全一致：
//   - 文件头 7 字节：":SWL" + uint16 小端魔数 0xE1D6 + uint8 版本号 2
//   - 其后为 LevelDB log 分帧：整个字节流按 32 KiB 分块，块内紧凑排列 chunk
//   - chunk 头 7 字节：CRC32(4B) + payload 长度(2B) + chunk 类型(1B)
//   - chunk 类型 1=FULL、2=FIRST、3=MIDDLE、4=LAST；跨块记录拆成 FIRST + MIDDLE* + LAST
//   - 块尾剩余不足 7 字节时以 0 填充，对齐到下一块
//   - CRC 校验范围是「类型字节 + payload」，即以 crc32(类型字节) 为预置基值
//
// payload 是序列化后的 swanlab.record.v1.Record。本包只负责字节搬运，
// 不 import proto，也不解释记录语义。
package store

import (
	"errors"
	"hash/crc32"
)

const (
	// HeaderLen 是文件头与 chunk 头的长度，两者恰好都是 7 字节。
	HeaderLen = 7
	// BlockLen 是物理块大小。
	BlockLen = 32768
	// HeaderVersion 是当前文件格式版本号，与 Python 侧保持一致。
	HeaderVersion = 2

	// dataLen 是单个 chunk 能承载的最大 payload 长度。
	dataLen = BlockLen - HeaderLen

	headerIdent = ":SWL"
	headerMagic = 0xE1D6

	chunkFull   byte = 1
	chunkFirst  byte = 2
	chunkMiddle byte = 3
	chunkLast   byte = 4
)

var (
	// ErrBadHeader 表示文件头不符合 :SWL 格式，或版本号不兼容。
	ErrBadHeader = errors.New("store: invalid file header")
	// ErrCorrupt 表示帧结构非法或 CRC 校验失败，其后的数据已不可信。
	ErrCorrupt = errors.New("store: corrupt record")
	// ErrTruncated 表示文件在记录中间结束（写入方尚未写完或进程异常退出）。
	ErrTruncated = errors.New("store: truncated record")
	// ErrClosed 表示操作发生在已关闭的 Writer 或 AsyncWriter 上。
	ErrClosed = errors.New("store: closed")
)

// crcSeed[t] 是类型字节 t 的 CRC32，作为该 chunk 校验和的预置基值。
var crcSeed = func() (s [5]uint32) {
	for t := chunkFull; t <= chunkLast; t++ {
		s[t] = crc32.ChecksumIEEE([]byte{t})
	}
	return s
}()
