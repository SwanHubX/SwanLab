// Code generated by protoc-gen-go. DO NOT EDIT.
// versions:
// 	protoc-gen-go v1.36.6
// 	protoc        v6.31.0
// source: swanlab/proto/core/collector/v1/collector_service.proto

package pb

import (
	protoreflect "google.golang.org/protobuf/reflect/protoreflect"
	protoimpl "google.golang.org/protobuf/runtime/protoimpl"
	reflect "reflect"
	sync "sync"
	unsafe "unsafe"
)

const (
	// Verify that this generated code is sufficiently up-to-date.
	_ = protoimpl.EnforceVersion(20 - protoimpl.MinVersion)
	// Verify that runtime/protoimpl is sufficiently up-to-date.
	_ = protoimpl.EnforceVersion(protoimpl.MaxVersion - 20)
)

type CollectorUploadRequest struct {
	state         protoimpl.MessageState `protogen:"open.v1"`
	Data          *KeyValueList          `protobuf:"bytes,1,opt,name=data,proto3" json:"data,omitempty"`
	unknownFields protoimpl.UnknownFields
	sizeCache     protoimpl.SizeCache
}

func (x *CollectorUploadRequest) Reset() {
	*x = CollectorUploadRequest{}
	mi := &file_swanlab_proto_core_collector_v1_collector_service_proto_msgTypes[0]
	ms := protoimpl.X.MessageStateOf(protoimpl.Pointer(x))
	ms.StoreMessageInfo(mi)
}

func (x *CollectorUploadRequest) String() string {
	return protoimpl.X.MessageStringOf(x)
}

func (*CollectorUploadRequest) ProtoMessage() {}

func (x *CollectorUploadRequest) ProtoReflect() protoreflect.Message {
	mi := &file_swanlab_proto_core_collector_v1_collector_service_proto_msgTypes[0]
	if x != nil {
		ms := protoimpl.X.MessageStateOf(protoimpl.Pointer(x))
		if ms.LoadMessageInfo() == nil {
			ms.StoreMessageInfo(mi)
		}
		return ms
	}
	return mi.MessageOf(x)
}

// Deprecated: Use CollectorUploadRequest.ProtoReflect.Descriptor instead.
func (*CollectorUploadRequest) Descriptor() ([]byte, []int) {
	return file_swanlab_proto_core_collector_v1_collector_service_proto_rawDescGZIP(), []int{0}
}

func (x *CollectorUploadRequest) GetData() *KeyValueList {
	if x != nil {
		return x.Data
	}
	return nil
}

type CollectorUploadResponse struct {
	state         protoimpl.MessageState `protogen:"open.v1"`
	Success       bool                   `protobuf:"varint,1,opt,name=success,proto3" json:"success,omitempty"`
	Message       string                 `protobuf:"bytes,2,opt,name=message,proto3" json:"message,omitempty"`
	unknownFields protoimpl.UnknownFields
	sizeCache     protoimpl.SizeCache
}

func (x *CollectorUploadResponse) Reset() {
	*x = CollectorUploadResponse{}
	mi := &file_swanlab_proto_core_collector_v1_collector_service_proto_msgTypes[1]
	ms := protoimpl.X.MessageStateOf(protoimpl.Pointer(x))
	ms.StoreMessageInfo(mi)
}

func (x *CollectorUploadResponse) String() string {
	return protoimpl.X.MessageStringOf(x)
}

func (*CollectorUploadResponse) ProtoMessage() {}

func (x *CollectorUploadResponse) ProtoReflect() protoreflect.Message {
	mi := &file_swanlab_proto_core_collector_v1_collector_service_proto_msgTypes[1]
	if x != nil {
		ms := protoimpl.X.MessageStateOf(protoimpl.Pointer(x))
		if ms.LoadMessageInfo() == nil {
			ms.StoreMessageInfo(mi)
		}
		return ms
	}
	return mi.MessageOf(x)
}

// Deprecated: Use CollectorUploadResponse.ProtoReflect.Descriptor instead.
func (*CollectorUploadResponse) Descriptor() ([]byte, []int) {
	return file_swanlab_proto_core_collector_v1_collector_service_proto_rawDescGZIP(), []int{1}
}

func (x *CollectorUploadResponse) GetSuccess() bool {
	if x != nil {
		return x.Success
	}
	return false
}

func (x *CollectorUploadResponse) GetMessage() string {
	if x != nil {
		return x.Message
	}
	return ""
}

var File_swanlab_proto_core_collector_v1_collector_service_proto protoreflect.FileDescriptor

const file_swanlab_proto_core_collector_v1_collector_service_proto_rawDesc = "" +
	"\n" +
	"7swanlab/proto/core/collector/v1/collector_service.proto\x12\x1fswanlab.proto.core.collector.v1\x1a$swanlab/proto/common/v1/common.proto\"S\n" +
	"\x16CollectorUploadRequest\x129\n" +
	"\x04data\x18\x01 \x01(\v2%.swanlab.proto.common.v1.KeyValueListR\x04data\"M\n" +
	"\x17CollectorUploadResponse\x12\x18\n" +
	"\asuccess\x18\x01 \x01(\bR\asuccess\x12\x18\n" +
	"\amessage\x18\x02 \x01(\tR\amessage2\x88\x01\n" +
	"\tCollector\x12{\n" +
	"\x06Upload\x127.swanlab.proto.core.collector.v1.CollectorUploadRequest\x1a8.swanlab.proto.core.collector.v1.CollectorUploadResponseB\rZ\vcore/pkg/pbb\x06proto3"

var (
	file_swanlab_proto_core_collector_v1_collector_service_proto_rawDescOnce sync.Once
	file_swanlab_proto_core_collector_v1_collector_service_proto_rawDescData []byte
)

func file_swanlab_proto_core_collector_v1_collector_service_proto_rawDescGZIP() []byte {
	file_swanlab_proto_core_collector_v1_collector_service_proto_rawDescOnce.Do(func() {
		file_swanlab_proto_core_collector_v1_collector_service_proto_rawDescData = protoimpl.X.CompressGZIP(unsafe.Slice(unsafe.StringData(file_swanlab_proto_core_collector_v1_collector_service_proto_rawDesc), len(file_swanlab_proto_core_collector_v1_collector_service_proto_rawDesc)))
	})
	return file_swanlab_proto_core_collector_v1_collector_service_proto_rawDescData
}

var file_swanlab_proto_core_collector_v1_collector_service_proto_msgTypes = make([]protoimpl.MessageInfo, 2)
var file_swanlab_proto_core_collector_v1_collector_service_proto_goTypes = []any{
	(*CollectorUploadRequest)(nil),  // 0: swanlab.proto.core.collector.v1.CollectorUploadRequest
	(*CollectorUploadResponse)(nil), // 1: swanlab.proto.core.collector.v1.CollectorUploadResponse
	(*KeyValueList)(nil),            // 2: swanlab.proto.common.v1.KeyValueList
}
var file_swanlab_proto_core_collector_v1_collector_service_proto_depIdxs = []int32{
	2, // 0: swanlab.proto.core.collector.v1.CollectorUploadRequest.data:type_name -> swanlab.proto.common.v1.KeyValueList
	0, // 1: swanlab.proto.core.collector.v1.Collector.Upload:input_type -> swanlab.proto.core.collector.v1.CollectorUploadRequest
	1, // 2: swanlab.proto.core.collector.v1.Collector.Upload:output_type -> swanlab.proto.core.collector.v1.CollectorUploadResponse
	2, // [2:3] is the sub-list for method output_type
	1, // [1:2] is the sub-list for method input_type
	1, // [1:1] is the sub-list for extension type_name
	1, // [1:1] is the sub-list for extension extendee
	0, // [0:1] is the sub-list for field type_name
}

func init() { file_swanlab_proto_core_collector_v1_collector_service_proto_init() }
func file_swanlab_proto_core_collector_v1_collector_service_proto_init() {
	if File_swanlab_proto_core_collector_v1_collector_service_proto != nil {
		return
	}
	file_swanlab_proto_common_v1_common_proto_init()
	type x struct{}
	out := protoimpl.TypeBuilder{
		File: protoimpl.DescBuilder{
			GoPackagePath: reflect.TypeOf(x{}).PkgPath(),
			RawDescriptor: unsafe.Slice(unsafe.StringData(file_swanlab_proto_core_collector_v1_collector_service_proto_rawDesc), len(file_swanlab_proto_core_collector_v1_collector_service_proto_rawDesc)),
			NumEnums:      0,
			NumMessages:   2,
			NumExtensions: 0,
			NumServices:   1,
		},
		GoTypes:           file_swanlab_proto_core_collector_v1_collector_service_proto_goTypes,
		DependencyIndexes: file_swanlab_proto_core_collector_v1_collector_service_proto_depIdxs,
		MessageInfos:      file_swanlab_proto_core_collector_v1_collector_service_proto_msgTypes,
	}.Build()
	File_swanlab_proto_core_collector_v1_collector_service_proto = out.File
	file_swanlab_proto_core_collector_v1_collector_service_proto_goTypes = nil
	file_swanlab_proto_core_collector_v1_collector_service_proto_depIdxs = nil
}
