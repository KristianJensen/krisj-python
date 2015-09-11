

from __future__ import absolute_import, print_function

import math
import bokeh.plotting as bplt
from bokeh.charts import Bar, Histogram, HeatMap
import colorsys

import numpy as np
import pandas as pd


def color_range(n, lightness=0.5, saturation=1, offset=0):
    hls_list = [(i+offset, lightness, saturation) for i in np.linspace(0, 1, n, endpoint=False)]
    rgb_list = [colorsys.hls_to_rgb(*ele) for ele in hls_list]
    hex_list = [generate_hex_string(*ele) for ele in rgb_list]
    return hex_list


def generate_hex_string(r, g, b):
    return "#%2.2x%2.2x%2.2x" % (r*255, g*255, b*255)


def dataframe_summary_plot(df, categorical_column, numerical_column, error_bars="std", **kwargs):
    categories = list(set(df[categorical_column]))
    n_cats = len(categories)
    xs = np.arange(1, n_cats+1)

    width = 0.3

    means = []
    stds = []
    sems = []
    for i, category in enumerate(categories):
        dat = df[df[categorical_column] == category][numerical_column]
        mean = dat.mean()
        std = dat.std(ddof=1)
        sem = std / math.sqrt(len(dat))

        means.append(mean)
        stds.append(std)
        sems.append(sem)

    means = np.array(means)
    stds = np.array(stds)
    sems = np.array(sems)


    if error_bars == "std":
        errors = stds
    elif error_bars == "sem":
        errors = sems
    else:
        errors = None

    lowers = means - errors
    uppers = means + errors


    data = {numerical_column: (means, lowers, uppers)}

    fig = bar_plot_with_error_bars(data, categories=categories, **kwargs)

    return fig


def bar_plot_with_error_bars(data, categories=None, colors=None, legend=True):
    """

    :param data: Dict containing series. Each series containing a tuple of (means, lowers, uppers) each of length n_cats.
    :return: bokeh figure
    """
    n_series = len(data)
    n_cats = None

    width = 0.5

    if colors is None:
        colors = color_range(n_series)

    for series, (means, lowers, uppers) in data.items():
        if len(means) == len(lowers) == len(uppers):
            if n_cats is None:
                n_cats = len(means)
            else:
                if len(means) != n_cats:
                    raise ValueError("Unequal number of categories in series")
        else:
            raise ValueError("Unequal number of categories in series")

    if categories is None:
        categories = list(map(str, range(1, n_cats+1)))
    else:
        if len(categories) != n_cats:
            raise ValueError("Length of 'categories' argument does not fit the data")

    fig = bplt.figure(x_range=categories)
    xs = np.arange(1, n_cats+1)

    series_centers = np.linspace(-0.5, 0.5, n_series+2)[1:-1]

    for i, series in enumerate(data):

        means, lowers, uppers = map(np.array, (data[series]))

        if legend:
            series_legend=series
        else:
            series_legend=None

        fig.quad(bottom=[0]*n_cats, top=means,
                 left=xs+series_centers[i]-width/(2*n_series),
                 right=xs+series_centers[i]+width/(2*n_series), color=colors[i], legend=series_legend)

        for j in range(n_cats):
            fig.line( [xs[j]+series_centers[i]]*2, [lowers[j], uppers[j]], color="black")

    fig.xgrid.grid_line_color = None
    fig.xaxis.major_label_orientation = 1

    return fig




