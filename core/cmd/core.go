package main

// @Title        core.go
// @Description  Package main implements a server for gRPC service.
// @Create       kaikai 2025/6/11 11:24

import (
	"flag"
	"fmt"
	"log"
	"net"

	"google.golang.org/grpc"

	"github.com/SwanHubX/SwanLab/core/internal/service"
	"github.com/SwanHubX/SwanLab/core/pkg/pb"
)

var (
	port = flag.Int("port", 0, "The server port")
)

func main() {
	flag.Parse()
	lis, err := net.Listen("tcp", fmt.Sprintf(":%d", *port))
	if err != nil {
		log.Fatalf("failed to listen: %v", err)
	}
	s := grpc.NewServer()
	pb.RegisterHealthServiceServer(s, &service.HealthService{})
	log.Printf("server listening at %v", lis.Addr())
	if err := s.Serve(lis); err != nil {
		log.Fatalf("failed to serve: %v", err)
	}
}
