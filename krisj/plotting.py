

from __future__ import absolute_import, print_function

import math
import bokeh.plotting as bplt
import colorsys

import numpy as np
import random


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


def bar_plot_with_error_bars(data, categories=None, colors=None, legend=True, **kwargs):
    """

    :param data: Dict containing series. Each series containing a tuple of (means, lowers, uppers) each of length n_cats.
    :return: bokeh figure
    """
    n_series = len(data)
    n_cats = None

    column_width = 0.4
    width = 0.4

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

    fig = bplt.figure(x_range=categories, **kwargs)
    xs = np.arange(1, n_cats+1)

    series_centers = np.linspace(-width, width, n_series+2)[1:-1]

    for i, series in enumerate(data):

        means, lowers, uppers = map(np.array, (data[series]))

        if legend:
            series_legend = series
        else:
            series_legend = None

        fig.quad(bottom=[0]*n_cats, top=means,
                 left=xs+series_centers[i]-column_width/(2*n_series),
                 right=xs+series_centers[i]+column_width/(2*n_series), color=colors[i], legend=series_legend)

        for j in range(n_cats):
            fig.line( [xs[j]+series_centers[i]]*2, [lowers[j], uppers[j]], color="black")

    fig.xgrid.grid_line_color = None
    fig.xaxis.major_label_orientation = 1

    return fig


def make_categorical_scatter(data, color="#669900", percentiles=(), **kwargs):
    """
    :param data:
    :return:
    """
    circle_kwargs = {}
    for k, v in kwargs.items():
        if k.startswith("circle_"):
            circle_kwargs[k[7:]] = v
            del kwargs[k]

    categories = list(data.keys())

    fig = bplt.figure(title=None, x_range=categories, **kwargs)

    fig.xaxis.major_label_text_font_size = "16px"

    for i, name in enumerate(categories):
        dat_array = np.array(data[name])
        if 50 in percentiles or 50. in percentiles:
            median = np.percentile(dat_array, 50)
            fig.line([i+0.8, i+1.2],[median, median], color="black")
        dat_list = list(data[name])
        x_list = [i+1+random.gauss(0, 0.1) for _ in dat_list]
        fig.circle(x_list, dat_list, fill_alpha=0.4, line_alpha=0,  fill_color=color, **circle_kwargs)

        for perc in percentiles:
            percentile = np.percentile(dat_array, perc)
            fig.line([i+0.9, i+1.1], [percentile, percentile], color="black")

    return fig


def radar_plot(df, title=None):
    fig = bplt.figure(tools=["reset", "pan", "wheel_zoom", "save"], height=600,
                      width=600, y_range=[-1.6, 1.6], x_range=[-1.6, 1.6], title=title)

    fig.grid.grid_line_color = None
    fig.axis.axis_line_color = None
    fig.axis.major_tick_line_color = None
    fig.axis.minor_tick_line_color = None
    fig.axis.major_label_text_color = None
    fig.outline_line_color = None

    labels = df.columns
    n = len(labels)
    series = df.index
    data = df.values
    ticksize = 0.02
    ticknum = 5

    color_strings = color_range(len(data), 0.8)

    spokes = []

    max_scale = 0.5

    # Draw spokes
    for i in range(n):
        end = (math.cos(2.*math.pi/n*i), math.sin(2.*math.pi/n*i))
        spokes.append(end)
        text_offset = len(labels[i])/16.*int(end[0] < 0)
        fig.text(1.1*end[0]-text_offset, 1.1*end[1]-0.05, text=labels[i])
        fig.line([0, end[0]], [0, end[1]])
        # Draw ticks on spokes
        for j in range(1, ticknum+1):
            pt = end[0]/ticknum*j, end[1]/ticknum*j
            up = pt[0]+end[1]*ticksize, pt[1]-end[0]*ticksize
            down = pt[0]-end[1]*ticksize, pt[1]+end[0]*ticksize
            text_loc = pt[0]-end[1]*3*ticksize-0.05, pt[1]+end[0]*3*ticksize-0.05
            fig.line([up[0], down[0]], [up[1], down[1]])
            fig.text(*text_loc, text=str(max_scale/ticknum*j), text_font_size="12px")

    # Draw data
    for j, d in enumerate(data):
        xs = [d[i]*spokes[i][0]/max_scale for i in range(len(spokes))]
        ys = [d[i]*spokes[i][1]/max_scale for i in range(len(spokes))]
        fig.patch(x=xs, y=ys, color=color_strings[j], line_width=3, fill_alpha=0.1, legend=series[j])

    #fig.circle(x=0, y=0, radius=1, fill_alpha=0)

    return fig
