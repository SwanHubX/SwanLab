# 这里提供了 SwanLab 所有支持的图表类型

每个图表类型被抽象为一个组件,但是在外部通过本部分提供的[ChartsDashboard](./ChartsDashboard.vue)组件统一调用，所需的就是传入chart配置  
> 一个ChartsContainer组件代表一个命名空间
作为数据驱动的图表仪表盘，每一个图表都接受一个图表配置。与此同时，为了与业务逻辑解藕，一些图表显示行为、图表配置将由外部传入，这一些原则将与全局组件的设计保持一致。

## 使用图像仪表盘

图像仪表盘及[ChartsDashboard](./ChartsDashboard.vue)组件，将接收图表组织配置（groups）以及一些其他的配置。
如果需要使用图像仪表盘，你需要确保通过props传入了以下内容：

1. `groups`: 代表图表组织配置
2. `updateChartStatus`: 更新图表状态的方法
3. `updateNamespaceStatus`: 更新图表配置的方法
4. `media`: 媒体查询对象
5. `getColor`: 获取颜色的方法
6. `defaultColor`: 默认颜色

接下来逐个说明这些props。

### groups

代表图表组织配置，包含namespace等内容，具体可以看对接的api文档，这里不再赘述。

### updateChartStatus

更新图表状态（置顶还是隐藏或者恢复正常状态），接受图表对象和状态吗，例如：

```js
/**
 * 更新chart状态，置顶、隐藏或正常显示
 * @param { object } chart 图表对象
 * @param { int } status 状态码，0为正常，1为置顶，-1为隐藏
 * @returns { Promise } Promise对象，最终可返回更新后的图表组织结构
 */
export const updateChartStatus = async (chart, status) => {
  const { data } = await http.patch('/chart/' + chart.id + '/status', {
    status
  })
  return data
}
```

### updateNamespaceStatus

更新命名空间状态（开启还是关闭），接收命名空间对象和开启关闭参数，例如：

```js
/**
 * 更新命名空间展开状态
 * @param { boolean } opened 是否展开
 * @param { object } namespace 命名空间对象
 * @returns { Promise }
 */
export const updateNamespaceStatus = (opened, namespace) => {
  return http.patch('/namespace/' + namespace.id + '/opened', {
    experiment_id: namespace.experiment_id?.id,
    project_id: namespace.project_id?.id,
    opened
  })
}
```

### media

获取媒体文件，输入文件名、实验id和数据名称，从后端获取，例如:

```js
/**
 * 获取媒体文件，获取blob对象
 * 返回promise，如果成功，返回blob对象，否则返回错误信息
 */
export const media = {
  /**
   * 获取媒体文件，获取blob对象
   * 返回promise，如果成功，返回blob对象，否则返回错误信息
   * @param { string } data 即文件名，即数据中的data字段
   * @param { string } experiment_id 实验id
   * @param { string } tag 数据名称
   * @returns { Promise<Blob> }
   */
  get: (data, experiment_id, tag) => {
    return new Promise((resolve, reject) => {
      http
        .get('/media/' + data, { params: { tag, experiment_id }, responseType: 'blob' })
        .then((res) => {
          resolve(res)
        })
        .catch((err) => {
          // 如果blob的type为application/json，说明是后端返回的错误信息，需要解析
          if (err.data.type === 'application/json') {
            const reader = new FileReader()
            reader.readAsText(err.data)
            reader.onload = () => {
              reject(JSON.parse(reader.result))
            }
          } else reject(err)
        })
    })
  }
}
```

### defaultColor

图表默认颜色，为十六进制颜色值，例如：`#000000`，主要用于展示音频等图表的颜色。

### getColor

获取颜色的方法，接收一个source名称和当前source的index，返回一个颜色值。

### 订阅与取消订阅

图表仪表盘需要订阅数据以展示，因此需要传入订阅和取消订阅的方法。前者用于获取数据，后者用于避免不必要的数据请求。

对于订阅函数而言，图表内部将会传入以下几个参数来期望外部实现数据订阅：
| 属性 | 类型 | 说明 |
| --- | --- | --- |
| sources | array<string> | 即某个图表的sources字段 |
| cid | any | 某个图表的id |
| callback | function | 图表内部的回调函数，当数据请求完毕时应该执行此函数 |

对于回调函数而言，他应该接收以下几个参数：
| 属性 | 类型 | 说明 |
| --- | --- | --- |
| source | string | 告知图表是sources字段中的哪个源发生了更新 |
| data | object | 更新后的数据，详见下方的data |
| error | any | 错误信息 |

## 图表组件

所有的图表组件都将有统一的props输入：
| 属性 | 类型 | 说明 |
| --- | --- | --- |
| title | String | 图表标题 |
| chart | Object | 图表配置 |
| index | Number | 图表在当前命名空间中的排序 |

所有的图表也将有通用的api，通过defineExpose暴露给父组件：
| 方法 | 类型 | 说明 |
| --- | --- | --- |
| render | Function | 渲染图表 |
| change | Function | 重渲染图表 |
| zoom | Function | 放大图表 |
| smooth | Function | 图表平滑 |

其中，render、change和zoom将接受一个参数，即图表数据，其数据格式为：

```js
const data = {
    "数据名称": {
        list:[
            {
              data: "数据",
              step: "步骤"，
              create_time: "创建时间"
            },
            ...
        ],
        ... // 一些其他的配置字段
    }
    ...
}
```

这些数据由组件使用者提供，通过依赖注入的方式为内部组件提供订阅数据、更新数据的能力。

本模块不再遵守vue单向数据流规范，但是也不鼓励暴露数据，而应该暴露api，通过api来完成逻辑的处理。

大致的组件模版可以输入`vue-chart`，由代码片段生成。
