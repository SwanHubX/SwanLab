"""
@author: caddiesnew
@file: metrics.py
@time: 2025-07-27 19:22:45
@description: 基于 echarts 实现部分 metrics 图表
"""

import importlib.util

import numpy as np
import pyecharts.options as opts
from pyecharts.charts import HeatMap, Line


def _check_sklearn():
    """Check if scikit-learn is installed."""
    if importlib.util.find_spec("sklearn") is None:
        raise TypeError("scikit-learn is required for this function. Please install it via `pip install scikit-learn`.")


def roc_curve(ground_truth, predictions, title=None) -> Line:
    """Plot Receiver operating characteristic (ROC) with AUC

    Parameters
    ----------
    ground_truth : ArrayLike. array-like of shape (n_samples,)
        True binary labels.
    predictions : ArrayLike. array-like of shape (n_samples,)
        Target scores, can either be probability estimates of the positive
            class, confidence values, or non-thresholded measure of decisions
    title : str or bool, optional
        Title for the chart. If True, will use default title with AUC.
        If string, will use that as title. If None or False, no title.
    """
    _check_sklearn()
    from sklearn import metrics

    fpr, tpr, _ = metrics.roc_curve(ground_truth, predictions)
    roc_auc = metrics.auc(fpr, tpr)

    # Diagonostic line
    diag_arr = fpr.copy()

    line = (
        Line()
        .add_xaxis(xaxis_data=fpr.tolist())
        .add_yaxis(
            series_name="TPR",
            y_axis=tpr.tolist(),
            is_smooth=True,
            is_symbol_show=False,
            label_opts=opts.LabelOpts(is_show=False),
            linestyle_opts=opts.LineStyleOpts(width=2),
        )
        .add_yaxis(
            series_name="",
            y_axis=diag_arr.tolist(),
            is_symbol_show=False,
            label_opts=opts.LabelOpts(is_show=False),
            linestyle_opts=opts.LineStyleOpts(width=2, type_="dashed"),
        )
    )

    # 处理标题
    if title is True:
        title_opts = opts.TitleOpts(title=f"ROC Curve (AUC = {roc_auc:.2f})", pos_left="center")
    elif isinstance(title, str):
        title_opts = opts.TitleOpts(title=title, pos_left="center")
    else:
        title_opts = None

    # 设置全局选项
    global_opts = {
        "xaxis_opts": opts.AxisOpts(type_="value", name="FPR"),
        "yaxis_opts": opts.AxisOpts(type_="value", name="TPR"),
        "legend_opts": opts.LegendOpts(is_show=False),
    }

    if title_opts:
        global_opts["title_opts"] = title_opts

    line.set_global_opts(**global_opts)
    return line


def pr_curve(ground_truth, predictions, title=None) -> Line:
    """Plot Precision-Recall(PR) Curve and AUC.

    Parameters
    ----------
    ground_truth : ArrayLike. array-like of shape (n_samples,)
        True binary labels.
    predictions : ArrayLike. array-like of shape (n_samples,)
        Target scores, can either be probability estimates of the positive
            class, confidence values, or non-thresholded measure of decisions
    title : str or bool, optional
        Title for the chart. If True, will use default title with AUC.
        If string, will use that as title. If None or False, no title.
    """
    _check_sklearn()
    from sklearn import metrics

    precision, recall, _ = metrics.precision_recall_curve(ground_truth, predictions)
    pr_auc = metrics.auc(recall, precision)

    line = (
        Line()
        .add_xaxis(xaxis_data=recall.tolist())
        .add_yaxis(
            series_name="Precision",
            y_axis=precision.tolist(),
            is_smooth=True,
            is_symbol_show=False,
            label_opts=opts.LabelOpts(is_show=False),
            linestyle_opts=opts.LineStyleOpts(width=2),
        )
    )

    # 处理标题
    if title is True:
        title_opts = opts.TitleOpts(title=f"PR Curve (AUC = {pr_auc:.2f})", pos_left="center")
    elif isinstance(title, str):
        title_opts = opts.TitleOpts(title=title, pos_left="center")
    else:
        title_opts = None

    # 设置全局选项
    global_opts = {
        "xaxis_opts": opts.AxisOpts(type_="value", name="Recall"),
        "yaxis_opts": opts.AxisOpts(type_="value", name="Precision"),
        "legend_opts": opts.LegendOpts(is_show=False),
    }

    if title_opts:
        global_opts["title_opts"] = title_opts

    line.set_global_opts(**global_opts)
    return line


def confusion_matrix(ground_truth, predictions, class_names, title=None) -> HeatMap:
    """Compute and Plot Confusion Matrix

    Parameters
    ----------
        ground_truth: ArrayLike: Ground truth (correct) target values.
        predictions: ArrayLike: Estimated targets as returned by a classifier.
        class_names: (list[str]): class names
        title : str or bool, optional
            Title for the chart. If True, will use default title.
            If string, will use that as title. If None or False, no title.
    """
    _check_sklearn()
    from sklearn import metrics

    cm = metrics.confusion_matrix(ground_truth, predictions)
    # normalized
    cm = cm.astype("float") / cm.sum(axis=1)[:, np.newaxis]
    # 注意：pyecharts中[x, y, value]的x对应列(预测)，y对应行(真实)
    # pyecharts的y轴是从下往上，而混淆矩阵是从上往下，需要反转y坐标
    data = []
    for i in range(len(class_names)):
        for j in range(len(class_names)):
            data.append([class_names[i], class_names[j], cm[i][j]])

    heatmap = (
        HeatMap()
        .add_xaxis(class_names)
        .add_yaxis(
            series_name="Confusion Matrix",
            yaxis_data=class_names,  # type: ignore
            value=data,
            label_opts=opts.LabelOpts(is_show=True, position="inside"),
        )
    )

    # 处理标题
    if title is True:
        title_opts = opts.TitleOpts(title="Confusion Matrix", pos_left="center")
    elif isinstance(title, str):
        title_opts = opts.TitleOpts(title=title, pos_left="center")
    else:
        title_opts = None

    # 设置全局选项
    global_opts = {
        "xaxis_opts": opts.AxisOpts(
            name="Prediction",
            type_="category",
            splitline_opts=opts.SplitLineOpts(
                is_show=True,
                linestyle_opts=opts.LineStyleOpts(color="#fff", width=2),
            ),
        ),
        "yaxis_opts": opts.AxisOpts(
            name="Ground truth",
            type_="category",
            splitline_opts=opts.SplitLineOpts(
                is_show=True,
                linestyle_opts=opts.LineStyleOpts(color="#fff", width=2),
            ),
        ),
        "legend_opts": opts.LegendOpts(is_show=False),
        "visualmap_opts": opts.VisualMapOpts(
            pos_right="right",
            pos_top="center",
            orient="vertical",
            min_=0,
            max_=int(cm.max()),
            is_piecewise=False,
            # 使用 seaborn 默认的 colormap
            range_color=[
                "#f7fbff",
                "#deebf7",
                "#c6dbef",
                "#9ecae1",
                "#6baed6",
                "#4292c6",
                "#2171b5",
                "#08519c",
                "#08306b",
            ],
        ),
    }

    if title_opts:
        global_opts["title_opts"] = title_opts

    heatmap.set_global_opts(**global_opts)
    return heatmap
