package main

// @Title        core.go
// @Description  Package main implements a server for gRPC service.
// @Create       kaikai 2025/6/11 11:24

import (
	"context"
	"flag"
	"fmt"
	"log"
	"net"

	"github.com/SwanHubX/SwanLab/core/pb"
	"google.golang.org/grpc"
	"google.golang.org/grpc/metadata"
)

var (
	port = flag.Int("port", 50051, "The server port")
)

// server is used to implement health server.
type server struct {
	pb.UnimplementedHealthServiceServer
}

func (s *server) Check(ctx context.Context, in *pb.HealthCheckRequest) (*pb.HealthCheckResponse, error) {
	md, ok := metadata.FromIncomingContext(ctx)
	if !ok {
		return nil, fmt.Errorf("missing metadata")
	}
	log.Printf("metadata: %v", md.Get("version"))
	log.Printf("SDK Version: %s", in.Version)
	trailer := metadata.Pairs("trailer-key", "val")
	grpc.SetTrailer(ctx, trailer)
	return &pb.HealthCheckResponse{
		Status:  true,
		Message: "ok",
	}, nil
}

func main() {
	flag.Parse()
	lis, err := net.Listen("tcp", fmt.Sprintf(":%d", *port))
	if err != nil {
		log.Fatalf("failed to listen: %v", err)
	}
	s := grpc.NewServer()
	pb.RegisterHealthServiceServer(s, &server{})
	log.Printf("server listening at %v", lis.Addr())
	if err := s.Serve(lis); err != nil {
		log.Fatalf("failed to serve: %v", err)
	}
}
