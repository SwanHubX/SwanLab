// 前端图表部分以数组排序，每个元素为一个group即section
// section数组下标代表排序，下标越小排序越靠前
// section内有charts字段，为数组，内包含一个个charts对象
// 开源版对象类型以Open为前缀

/**
 * @typedef {Object} OpenSection 一个节容器，存放图表，和一些默认情况
 * @property {number | string} id 平常不会使用，用于区别pin的section（-1）和hidden的section（-2），不影响排序，只影响名称显示；另一方面，将使用这个作为v-for时的key
 * @property {string} name section的名称，特别的，如果为default或者id为-1、-2，则按照前端定义的名称显示
 * @property {boolean | number} opened section是否展开，默认为true或者1
 * @property {OpenChart[]} charts 这个section内包含的图表
 */

/**
 * @typedef {Object} OpenChart 一个图表配置，用于标注这个图表是什么类型以及配置是什么，使用了什么数据源
 * @property {number | string} id 图表唯一id，用于v-for时的key
 * @property {boolean} multi 标注是否为多实验图表显示，这涉及到图表组件内部不同的数据处理方式
 * @property {string} name 图表名称
 * @property {string} reference 参考系，目前默认为step
 * @property {string[]} source 数据源列表，这里是显示在图表上的数据名称
 * @property {Object} source_map 数据源映射，用于前端显示数据源的名称映射到数据源对应的实验id
 * @property {Object} error 当前图表是否错误，如果不为null说明图表错误，此信息用于显示错误
 * @property {string} error.data_class sdk中上传时的数据类型，用于前端提示
 */

/**
 * @typedef {Object} OpenMetricData 一个指标包含的数据
 * @property {string | number} experiment_id 实验id
 * @property {OpenMetricDetail[]} list 数据列表
 */

/**
 * @typedef {Object} OpenMetricDetail 一个指标包含的数据
 * @property {number} index 指标步数
 * @property {number | string | string[]} data 这一步的数据
 * @property {boolean} _last 是否为最后一条数据，最后一个数据设置为true即可，其他不需要设置
 */

/**
 * @typedef {Object<string, OpenMetricData>} OpenChartData 传递给图表的数据
 */

/**
 * @callback OpenChartSubscribe 图表本身向父组件传递数据的订阅函数
 * @param {string[]} sources 数据源名称，代表这个图表使用了哪些数据源
 * @param {string} chartId 图表id
 * @param {OpenChartSubscribeCallback} callback 订阅成功后的回调函数
 */

/**
 * @callback OpenChartSubscribeCallback 图表本身向父组件传递数据的 订阅成功/数据更改 的回调函数，父组件将在订阅成功后调用这个函数
 * @param {string} key 数据源名称
 * @param {string} data 数据源的数据
 * @param {Object} error 请求失败的错误信息——如果有的话，没有就是null
 */

// ---------------------------------- 图表api ----------------------------------

/**
 * @callback setOriginalChartHeight 设置图表的原始高度，如果在执行前图表还没有被渲染，那么将在图表渲染后立即执行
 * @param {number} height 图表的高度, 单位px
 * @param {number} [maxHeight=800] 图表的最大高度, 单位px, 默认为800
 * @param {number} [minHeight=200] 图表的最小高度, 单位px, 默认为200
 */
