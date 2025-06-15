// @Title        collector.go
// @Description  Package service implements the Collector service for gRPC.
// @Create       kaikai 2025/6/14 18:00

package service

import (
	"log/slog"

	"context"

	"github.com/SwanHubX/SwanLab/core/pkg/pb"
)

type Collector struct {
	pb.UnimplementedCollectorServer
}

func (c *Collector) Upload(ctx context.Context, in *pb.CollectorUploadRequest) (*pb.CollectorUploadResponse, error) {
	slog.InfoContext(ctx, "collector: Upload: started", "data", in.GetData())
	return &pb.CollectorUploadResponse{
		Success: true,
		Message: "ok",
	}, nil
}
