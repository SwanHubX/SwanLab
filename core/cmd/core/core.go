// @Title        core.go
// @Description  Package main implements a server for gRPC service.
// @Create       kaikai 2025/6/11 11:24

package main

import (
	"flag"
	"fmt"
	"log/slog"
	"net"

	"google.golang.org/grpc"
	"google.golang.org/grpc/health"
	grpcHealth "google.golang.org/grpc/health/grpc_health_v1"

	"github.com/SwanHubX/SwanLab/core/internal/service"
	"github.com/SwanHubX/SwanLab/core/pkg/pb"
)

func main() {
	port := flag.Int("port", 0, "The server port")
	flag.Parse()

	lis, err := net.Listen("tcp", fmt.Sprintf(":%d", *port)) //nolint:noctx  未来使用context代替
	if err != nil {
		slog.Error("failed to listen", "error", err)
	}

	s := grpc.NewServer()
	// Health Check Service
	healthCheck := health.NewServer()
	grpcHealth.RegisterHealthServer(s, healthCheck)
	// Collector Service
	collector := &service.Collector{}
	pb.RegisterCollectorServer(s, collector)
	healthCheck.SetServingStatus("collector", grpcHealth.HealthCheckResponse_SERVING)
	slog.Info("server listening at", "address", lis.Addr())

	if err = s.Serve(lis); err != nil {
		slog.Error("failed to serve", "error", err)
	}
}
