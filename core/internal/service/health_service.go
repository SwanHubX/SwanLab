package service

import (
	"context"
	"fmt"
	"log"

	"github.com/SwanHubX/SwanLab/core/pkg/pb"
	"google.golang.org/grpc/metadata"
)

// healthService is used to implement HealthServiceServer.
type HealthService struct {
	pb.UnimplementedHealthServiceServer
}

func (s *HealthService) Check(ctx context.Context, in *pb.HealthCheckRequest) (*pb.HealthCheckResponse, error) {
	md, ok := metadata.FromIncomingContext(ctx)
	if !ok {
		return nil, fmt.Errorf("missing metadata")
	}
	log.Printf("metadata: %v", md.Get("version"))
	log.Printf("SDK Version: %s", in.Version)
	return &pb.HealthCheckResponse{
		Status:  true,
		Message: "ok",
	}, nil
}
