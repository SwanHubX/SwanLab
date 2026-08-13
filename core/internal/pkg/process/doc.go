// Package process 管理 SwanLab core 与其父进程之间的生命周期关系。
//
// SwanLab core 通常由 Python SDK 作为子进程启动。正常退出时应由上层通过业务接口
// 通知 core 完成清理；本包处理的是父进程崩溃、被强制终止等无法发送退出通知的兜底场景，
// 防止 core 成为孤儿进程。
//
// # 设计原则
//
// 调用方在 core 启动阶段接收 Python SDK 传入的父进程 PID，并调用 NotifyOnParentExit 注册监控：
//
//	parentExited, err := process.NotifyOnParentExit(parentPID)
//	if err != nil {
//		// 监控未建立，调用方应终止启动，避免留下无法回收的 core 进程。
//	}
//	go func() {
//		<-parentExited
//		// 调用方在这里执行 flush、gRPC GracefulStop 等退出流程。
//	}()
//
// 注册前会验证 parentPID，平台实现完成初始化后还会再次检查父进程，避免父进程恰好在
// “检查存活”和“建立监控”之间退出的竞态。初始化失败会返回错误并通过 console 记录；
// 初始化成功及检测到父进程退出会记录 debug 日志。
//
// 本包只负责检测父进程退出并关闭通知 channel，不直接终止当前进程，也不执行任何业务清理。
// 退出顺序、清理超时和最终是否强制退出均由调用方决定，使正常退出与父进程异常退出可以复用
// 同一套生命周期编排。
//
// # 平台实现
//
//   - Linux：通过 prctl(PR_SET_PDEATHSIG, SIGUSR1) 让内核在父进程退出时发送可捕获信号，
//     再将该信号转换为 channel 通知。该路径不依赖轮询。SIGUSR1 由本包监听；收到信号后
//     仍会检查父 PID，其他来源的同名信号不会触发父进程退出通知。
//   - Windows：以 SYNCHRONIZE 权限打开父进程句柄，并用 WaitForSingleObject 等待句柄变为
//     signaled。该方式不依赖 Windows 上不会随父进程退出而变化的 PPID。
//   - 其他平台：定期比较 os.Getppid() 与预期 PID；父进程退出并发生 reparent 后关闭通知 channel。
//
// 本包位于 internal 下，仅供 SwanLab core 内部使用。
package process
