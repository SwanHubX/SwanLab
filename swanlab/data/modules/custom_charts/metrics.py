"""
@author: caddiesnew
@file: metrics.py
@time: 2025-07-27 19:22:45
@description: 基于 echarts 实现部分 metrics 图表
"""

import numpy as np
import pyecharts.options as opts
from pyecharts.charts import HeatMap, Line
from sklearn import metrics


def roc_curve(ground_truth, predictions) -> Line:
    """Plot Receiver operating characteristic (ROC) with AUC

    Parameters
    ----------
    ground_truth : ArrayLike. array-like of shape (n_samples,)
        True binary labels.
    predictions : ArrayLike. array-like of shape (n_samples,)
        Target scores, can either be probability estimates of the positive
            class, confidence values, or non-thresholded measure of decisions
    """
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
        .set_global_opts(
            xaxis_opts=opts.AxisOpts(type_="value", name="FPR"),
            yaxis_opts=opts.AxisOpts(type_="value", name="TPR"),
            legend_opts=opts.LegendOpts(is_show=False),
            # title_opts=opts.TitleOpts(title=f"PR Curve (AUC = {roc_auc:.2f})", pos_left="center"),
        )
    )
    return line


def pr_curve(ground_truth, predictions) -> Line:
    """Plot Precision-Recall(PR) Curve and AUC.

    Parameters
    ----------
    ground_truth : ArrayLike. array-like of shape (n_samples,)
        True binary labels.
    predictions : ArrayLike. array-like of shape (n_samples,)
        Target scores, can either be probability estimates of the positive
            class, confidence values, or non-thresholded measure of decisions
    """
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
        .set_global_opts(
            xaxis_opts=opts.AxisOpts(type_="value", name="Recall"),
            yaxis_opts=opts.AxisOpts(type_="value", name="Precision"),
            legend_opts=opts.LegendOpts(is_show=False),
            # title_opts=opts.TitleOpts(title=f"PR Curve (AUC = {pr_auc:.2f})", pos_left="center"),
        )
    )
    return line


def confusion_matrix(ground_truth, predictions, class_names: list[str]) -> HeatMap:
    """Compute and Plot Confusion Matrix

    Parameters
    ----------
        ground_trueth: ArrayLike: Ground truth (correct) target values.
        predictions: ArrayLike: Estimated targets as returned by a classifier.
        class_names: (list[str]): class names

    """
    cm = metrics.confusion_matrix(ground_truth, predictions)
    # normalzied
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
        .set_global_opts(
            xaxis_opts=opts.AxisOpts(
                name="Prediction",
                type_="category",
                splitline_opts=opts.SplitLineOpts(
                    is_show=True,
                    linestyle_opts=opts.LineStyleOpts(color="#fff", width=2),
                ),
            ),
            yaxis_opts=opts.AxisOpts(
                name="Ground truth",
                type_="category",
                splitline_opts=opts.SplitLineOpts(
                    is_show=True,
                    linestyle_opts=opts.LineStyleOpts(color="#fff", width=2),
                ),
            ),
            legend_opts=opts.LegendOpts(is_show=False),
            visualmap_opts=opts.VisualMapOpts(
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
        )
    )
    return heatmap
