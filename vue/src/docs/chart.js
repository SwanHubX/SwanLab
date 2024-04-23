// 前端图表部分以数组排序，每个元素为一个group即section
// section数组下标代表排序，下标越小排序越靠前
// section内有charts字段，为数组，内包含一个个charts对象
// 开源版对象类型以Open为前缀

/**
 * @description 一个完整的图表配置，包含了当前页面的所有图表信息
 * @typedef {OpenSection[]} OpenSections
 */

/**
 * @description 一个节容器，存放图表，和一些默认情况
 * @typedef {Object} OpenSection
 * @property {number | string} id 平常不会使用，用于区别pin的section（-1）和hidden的section（-2），不影响排序，只影响名称显示；另一方面，将使用这个作为v-for时的key
 * @property {string} name section的名称，特别的，如果为default或者id为-1、-2，则按照前端定义的名称显示
 * @property {boolean | number} opened section是否展开，默认为true或者1
 * @property {OpenChart[]} charts 这个section内包含的图表
 */

/**
 * @description 一个图表配置，用于标注这个图表是什么类型以及配置是什么，使用了什么数据源
 * @typedef {Object} OpenChart
 * @property {number | string} id 图表唯一id，用于v-for时的key
 * @property {boolean} multi 标注是否为多实验图表显示，这涉及到图表组件内部不同的数据处理方式
 * @property {string} name 图表名称
 * @property {string} reference 参考系，目前默认为step
 * @property {string[]} source 数据源列表，这里是显示在图表上的数据名称
 * @property {Object} source_map 数据源映射，用于前端显示数据源的名称映射到数据源对应的实验id
 * @property {Object} error 当前图表是否错误，如果不为null说明图表错误，此信息用于显示错误
 * @property {string} error.data_class sdk中上传时的数据类型，用于前端提示
 */
