# SwanLab API 文档

本目录存放 SwanLab 开放接口的 OpenAPI/Swagger 2.0 定义，合并了原 SwanLab-Auth 与 SwanLab-Server 两个仓库的规范。

## 目录内容

| 文件         | 说明            |
| ------------ | --------------- |
| `zh-cn.json` | 中文版 API 定义 |
| `en.json`    | 英文版 API 定义 |

> ⚠️ 修改时需保持两个文件**结构完全一致**（路径、字段、枚举、错误码均相同），仅 `description`、`summary` 等描述性文字允许翻译差异。

## 接口地址

所有接口共用 `host: api.swanlab.cn`、`basePath: /api`，完整 URL 为 `https://api.swanlab.cn/api` 加上各自的路径：

| 路径前缀          | 完整示例                                            |
| ----------------- | --------------------------------------------------- |
| `/auth/oauth2/*`  | `https://api.swanlab.cn/api/auth/oauth2/token`      |
| `/self_hosted/*`  | `https://api.swanlab.cn/api/self_hosted/root`       |

其中 `/auth/oauth2/*` 对应 SwanLab-Auth 原有的 `basePath: /api/auth/oauth2`，因 Swagger 2.0 单文件只允许一个顶层 `basePath`，故并入 `basePath: /api` 后写成路径段；接口完整 URL 与线上一致，未做改动。

- `securityDefinitions` 保留 `basicAuth`、`bearerAuth`、`apiKeyAuth` 三套认证方案。
- 顶层 `tags` 分组：`OAuth2`（OAuth2 授权与使用）与 `Self-hosted`（私有化部署管理）。

## 校验流程

改动 `*.json` 后，按以下流程校验。

### 1. 规范校验（必做）

校验文件是否符合 Swagger 2.0 规范：

```bash
npx --yes @redocly/cli lint docs/swagger/zh-cn.json docs/swagger/en.json
```

### 2. 渲染预览（建议）

按 Swagger UI 外观预览单个文档：

```bash
docker run --rm -p 8080:8080 \
  -v "$PWD/docs/swagger:/specs" \
  -e SWAGGER_JSON=/specs/zh-cn.json \
  swaggerapi/swagger-ui
```

打开 `http://localhost:8080` 查看；切换英文版把 `SWAGGER_JSON` 改成 `/specs/en.json` 即可。
