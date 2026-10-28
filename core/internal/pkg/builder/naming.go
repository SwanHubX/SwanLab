package builder

import "strconv"

// 本文件提供默认实验名与颜色词表：按历史实验数取模，
// 词表内容与顺序必须与 Python 客户端一致。

// presetColors 是 16 个预设颜色，下标即取模基。
var presetColors = [16]string{
	"#528d59", // 绿色
	"#587ad2", // 蓝色
	"#c24d46", // 红色
	"#9cbe5d", // 青绿色
	"#6ebad3", // 天蓝色
	"#dfb142", // 橙色
	"#6d4ba4", // 紫色
	"#8cc5b7", // 淡青绿色
	"#892d58", // 紫红色
	"#40877c", // 深青绿色
	"#d0703c", // 深橙色
	"#d47694", // 粉红色
	"#e3b292", // 淡橙色
	"#b15fbb", // 浅紫红色
	"#905f4a", // 棕色
	"#989fa3", // 灰色
}

// presetAnimals 是 20 个动物名词。
var presetAnimals = [20]string{
	"swan", "rat", "ox", "tiger", "rabbit",
	"dragon", "snake", "horse", "goat", "monkey",
	"rooster", "dog", "pig", "cat", "elephant",
	"penguin", "kangaroo", "panda", "lion", "zebra",
}

// normalAdjectives 是 20 个普通形容词。
var normalAdjectives = [20]string{
	"quick", "slow", "small", "large", "young",
	"old", "quiet", "loud", "happy", "calm",
	"busy", "lazy", "clever", "simple", "steady",
	"brave", "shy", "warm", "cool", "fresh",
}

// GenerateColor 按序号取模返回预设颜色：PRESET_COLORS[slug % 16]。
func GenerateColor(slug int) string {
	return presetColors[pyMod(slug, len(presetColors))]
}

// GenerateName 按序号取模生成默认实验名：
//
//	NORMAL_ADJECTIVES[slug % 20]-{PRESET_ANIMALS[slug % 20]}-{slug}
//
// 例如 slug=5 → "old-dragon-5"、slug=128 → "happy-goat-128"（与 Python 一致）。
func GenerateName(slug int) string {
	i := pyMod(slug, len(presetAnimals))
	return normalAdjectives[i] + "-" + presetAnimals[i] + "-" + strconv.Itoa(slug)
}
