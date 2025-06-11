package service

import (
	"log/slog"

	"golang.org/x/net/context"

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
