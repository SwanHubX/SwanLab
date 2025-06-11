package main

// @Title        core.go
// @Description  Package main implements a server for gRPC service.
// @Create       kaikai 2025/6/11 11:24

import (
	"flag"
	"fmt"
	"log/slog"
	"net"

	"google.golang.org/grpc"

	"github.com/SwanHubX/SwanLab/core/internal/service"
	"github.com/SwanHubX/SwanLab/core/pkg/pb"
)

func main() {
	port := flag.Int("port", 0, "The server port")
	flag.Parse()

	lis, err := net.Listen("tcp", fmt.Sprintf(":%d", *port))
	if err != nil {
		slog.Error("failed to listen", "error", err)
	}

	s := grpc.NewServer()
	pb.RegisterHealthServiceServer(s, &service.HealthService{})
	slog.Info("server listening at", "address", lis.Addr())

	if err = s.Serve(lis); err != nil {
		slog.Error("failed to serve", "error", err)
	}
}
