package service

import (
	"context"
	"errors"
	"log/slog"

	"google.golang.org/grpc/metadata"

	"github.com/SwanHubX/SwanLab/core/pkg/pb"
)

// HealthService is used to implement HealthServiceServer.
type HealthService struct {
	pb.UnimplementedHealthServiceServer
}

func (s *HealthService) Check(ctx context.Context, in *pb.HealthCheckRequest) (*pb.HealthCheckResponse, error) {
	md, ok := metadata.FromIncomingContext(ctx)
	if !ok {
		return nil, errors.New("missing metadata")
	}
	slog.InfoContext(ctx, "metadata", "version", md.Get("version"))
	slog.InfoContext(ctx, "request", "version", in.GetVersion())
	return &pb.HealthCheckResponse{
		Status:  true,
		Message: "ok",
	}, nil
}
