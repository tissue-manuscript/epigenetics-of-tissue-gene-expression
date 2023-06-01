#!/usr/bin/env python3
import csv
import math
import subprocess
import plotly.graph_objects as go
from gooey import Gooey, GooeyParser


def round_to_nearest(number: float, ref: float):
    return math.ceil(number / ref) * ref


def get_delimiter(input: str):
    if input == "comma":
        return ","
    elif input == "tab":
        return "\t"
    elif input == "commaspace":
        return ", "
    else:
        raise ValueError("Invalid delimiter specified!")


def produce_dot_plot(input_csv: str, input_delimiter: str, min_marker_size: int, max_marker_size: int,
                     interest_color_rgb: list, non_interest_color_rgb: list, grid_color_rgb: list,
                     grid_opacity: float, figure_width: int, figure_height: int, output_file: str, output_scale: float,
                     axis_tick_color_rgb: list, x_axis_start: float, x_axis_end: float, x_axis_tick_interval: int,
                     x_axis_title: str, x_axis_title_font_size: int, x_axis_tick_font_size: int, y_axis_tick_font_size: int):
    pvalue_log_function = lambda pv: -math.log10(pv)

    marker_size_range = max_marker_size - min_marker_size
    all_info = []

    with open(input_csv, "r") as data_file:
        reader = csv.reader(data_file, delimiter=input_delimiter)
        next(reader)  # skip header
        for row in reader:
            d = {}
            d["pathway"] = row[0]
            d["lpv"] = pvalue_log_function(float(row[1]))
            d["i"] = bool(int(row[2]))
            all_info.append(d)

    all_lpvs = [d["lpv"] for d in all_info]
    min_lpv, max_lpv = min(all_lpvs), max(all_lpvs)
    lpv_range = max_lpv - min_lpv

    new_all_info = sorted(all_info, key=lambda d: d["lpv"])

    fig = go.Figure()

    for d in new_all_info:
        print(d)

        marker_size = (((d["lpv"] - min_lpv) / lpv_range) * marker_size_range) + min_marker_size
        interest_color = f'rgba({interest_color_rgb[0]}, {interest_color_rgb[1]}, {interest_color_rgb[2]}, 0.95)'
        non_interest_color = \
            f'rgba({non_interest_color_rgb[0]}, {non_interest_color_rgb[1]}, {non_interest_color_rgb[2]}, 0.95)'

        fig.add_trace(go.Scatter(
            mode='markers',
            x=[d["lpv"]],
            y=[d["pathway"]],
            marker=dict(
                color=interest_color if d["i"] else non_interest_color,
                line_color='rgba(156, 165, 196, 1.0)',  # TODO
                line_width=1,
                symbol='circle',
                size=marker_size
            )
        ))

    grid_color = f'rgba({grid_color_rgb[0]}, {grid_color_rgb[1]}, {grid_color_rgb[2]}, {grid_opacity})'
    axis_tick_color = f'rgba({axis_tick_color_rgb[0]}, {axis_tick_color_rgb[1]}, {axis_tick_color_rgb[2]})'

    fig.update_layout(
        xaxis=dict(
            showgrid=True,
            showline=True,
            linecolor=axis_tick_color,
            tickfont_color=axis_tick_color,
            showticklabels=True,
            range=[x_axis_start, round_to_nearest(max_lpv, 0.5) if x_axis_end is None else x_axis_end],
            ticks='inside',
            dtick=x_axis_tick_interval,
            tickcolor=axis_tick_color,
            mirror=True,
            gridcolor=grid_color
        ),
        yaxis=dict(
            showgrid=True,
            showline=True,
            linecolor=axis_tick_color,
            tickfont_color=axis_tick_color,
            showticklabels=True,
            ticks='outside',
            tickcolor=axis_tick_color,
            mirror=True,
            gridcolor=grid_color
        ),
        margin=dict(l=140, r=40, b=50, t=80),
        showlegend=False,
        width=figure_width,
        height=int(len(all_info) / 5) * 300 if figure_height is None else figure_height,
        paper_bgcolor='white',
        plot_bgcolor='white',
        hovermode='closest',
    )

    if x_axis_title is not None and len(x_axis_title) > 0:
        fig.update_layout(
            xaxis_title=x_axis_title
        )
        fig.update_xaxes(
            title_font=dict(
                size=x_axis_title_font_size,
                color="Black"
            ))

    fig.update_xaxes(
        tickfont=dict(
            size=x_axis_tick_font_size
        )
    )

    fig.update_yaxes(
        tickfont=dict(
            size=y_axis_tick_font_size
        )
    )

    print("Writing figure to: " + output_file)
    fig.write_image(output_file, scale=output_scale)
    print ("Done!")
    # subprocess.run(['xdg-open', output_file], check=True)


@Gooey
def main():
    parser = GooeyParser(description='Produce a dot plot.')

    input_group = parser.add_argument_group(
        "Input File",
        "Input file options."
    )
    input_group.add_argument("-i", "--input", action="store", default="sample.csv", metavar="INPUT_CSV",
                             dest="input_file", widget='FileChooser',
                             help="Path to the input CSV (by default, 'sample.csv' in the current directory).")
    input_group.add_argument("--csvdelim", action="store", default="comma", metavar="CSV_DELIMITER",
                             choices=["comma", "tab", "commaspace"],
                             dest="file_delimiter",
                             help="The delimiter for the input file (default: 'comma', also valid: 'tab', 'commaspace')")

    output_group = parser.add_argument_group(
        "Output File",
        "Output file options."
    )
    output_group.add_argument("-o", "--output", action="store", default="fig.png", metavar="OUTPUT_PATH",
                              dest="output_file", widget='FileSaver',
                              help="Path to the output png (by default, 'fig.png' in the current directory).")
    output_group.add_argument("--scale", action="store", default=2, metavar="OUTPUT_FIG_SCALE",
                              dest="output_scale", type=float,
                              help="How much to scale the resolution of the output figure (default: 2).")

    graph_group = parser.add_argument_group(
        "Graph",
        "Graph options."
    )
    graph_group.add_argument("--min_marker_size", action="store", default=10, metavar="MIN_MARKER_SIZE", type=int,
                             dest="min_marker_size", help="The minimum size for a dot marker (default: 10)")
    graph_group.add_argument("--max_marker_size", action="store", default=40, metavar="MAX_MARKER_SIZE", type=int,
                             dest="max_marker_size", help="The maximum size for a dot marker (default: 40)")
    graph_group.add_argument("--interest-color-rgb", action="store", nargs=3, default='153 0 0', type=int,
                             dest="interest_color_rgb",
                             help="The RGB values for a dot of interest [default: 153 0 0]. " +
                                  "Valid values: Only integers between 0 and 255, inclusive.",
                             gooey_options={
                                 'validator': {
                                     'test': 'all(int(i) <= 255 and int(i) >= 0 for i in user_input.split(" ")) and len(user_input.split(" ")) == 3',
                                     'message': 'Must have 3 elements. Each element must be between 0 and 255, inclusive.'
                                 }
                             }
                             )
    graph_group.add_argument("--non-interest-color-rgb", action="store", nargs=3, default='0 0 179', type=int,
                             dest="non_interest_color_rgb",
                             help="The RGB values for a dot of interest [default: 0 0 179]. " +
                                  "Valid values: Only integers between 0 and 255, inclusive.",
                             gooey_options={
                                 'validator': {
                                     'test': 'all(int(i) <= 255 and int(i) >= 0 for i in user_input.split(" ")) and len(user_input.split(" ")) == 3',
                                     'message': 'Must have 3 elements. Each element must be between 0 and 255, inclusive.'
                                 }
                             }
                             )
    graph_group.add_argument("--grid-color-rgb", action="store", nargs=3, default='220 220 200', type=int,
                             dest="grid_color_rgb",
                             help="The RGB values for the background grid [default: 220 220 200]. " +
                                  "Valid values: Only integers between 0 and 255, inclusive.",
                             gooey_options={
                                 'validator': {
                                     'test': 'all(int(i) <= 255 and int(i) >= 0 for i in user_input.split(" ")) and len(user_input.split(" ")) == 3',
                                     'message': 'Must have 3 elements. Each element must be between 0 and 255, inclusive.'
                                 }
                             }
    )
    graph_group.add_argument("--axis-tick-color-rgb", action="store", nargs=3, default='0 0 0', type=int,
                             dest="axis_tick_color_rgb",
                             help="The RGB values for the axes and ticks [default: 0 0 0]. " +
                                  "Valid values: Only integers between 0 and 255, inclusive.",
                             gooey_options={
                                 'validator': {
                                     'test': 'all(int(i) <= 255 and int(i) >= 0 for i in user_input.split(" ")) and len(user_input.split(" ")) == 3',
                                     'message': 'Must have 3 elements. Each element must be between 0 and 255, inclusive.'
                                 }
                             }
                             )
    graph_group.add_argument("--grid_color_opacity", action="store", default=0.35, metavar="GRID_OPACITY", type=float,
                             dest="grid_opacity", help="The opacity of the background grid (default: 0.35). " +
                                                       "Valid values: Only floats between 0 and 1, inclusive",
                             gooey_options={
                                 'validator': {
                                     'test': 'float(user_input) >= 0.0 and float(user_input) <= 1.0',
                                     'message': 'Opacity must be between 0 and 1, inclusive.'
                                 }
                             }
                             )
    graph_group.add_argument("--x_axis_tick_interval", action="store", default=0.1, metavar="X_AXIS_TICK_INT", type=float,
                             dest="x_axis_tick_int", help="The interval between ticks on the x-axis (default: 0.1).")

    graph_group.add_argument("--width", action="store", default=500, metavar="FIGURE_WIDTH", type=int,
                             dest="fig_width", help="The width of the figure in pixels (default: 500).")
    graph_group.add_argument("--height", action="store", default=None, metavar="FIGURE_HEIGHT", type=int,
                             dest="fig_height", help="The height of the figure in pixels (default: blank = automatically calculated).")

    graph_group.add_argument("--x_axis_start", action="store", default=0.00, metavar="X_AXIS_START", type=float,
                             dest="x_axis_start", help="The starting point of the x-axis (default: 0.0).")
    graph_group.add_argument("--x_axis_ending", action="store", default=None, metavar="X_AXIS_END", type=float,
                             dest="x_axis_end", help="The ending point of the x-axis (default: blank = automatically calculated).")

    graph_group.add_argument("--x_axis_title", action="store", default=None, metavar="X_AXIS_TITLE",
                             type=str,
                             dest="x_axis_title", help="The label of the x-axis (default: blank).")
    graph_group.add_argument("--x_axis_title_font_size", action="store", default=12, metavar="X_AXIS_TITLE_FONT_SIZE",
                             type=int,
                             dest="x_axis_title_font_size", help="The font size of the x-axis label (default: 12).")
    graph_group.add_argument("--x_axis_tick_font_size", action="store", default=12, metavar="X_AXIS_TICK_FONT_SIZE",
                             type=int,
                             dest="x_axis_tick_font_size", help="The font size of the x-axis ticks (default: 12).")
    graph_group.add_argument("--y_axis_tick_font_size", action="store", default=12, metavar="Y_AXIS_TICK_FONT_SIZE",
                             type=int,
                             dest="y_axis_tick_font_size", help="The font size of the y-axis ticks (default: 12).")


    args = parser.parse_args()

    produce_dot_plot(input_csv=args.input_file, input_delimiter=get_delimiter(args.file_delimiter),
                     min_marker_size=args.min_marker_size, max_marker_size=args.max_marker_size,
                     interest_color_rgb=args.interest_color_rgb, non_interest_color_rgb=args.non_interest_color_rgb,
                     grid_color_rgb=args.grid_color_rgb, grid_opacity=args.grid_opacity,
                     figure_width=args.fig_width, figure_height=args.fig_height, output_file=args.output_file,
                     output_scale=args.output_scale, axis_tick_color_rgb=args.axis_tick_color_rgb,
                     x_axis_start=args.x_axis_start, x_axis_end=args.x_axis_end, x_axis_tick_interval=args.x_axis_tick_int,
                     x_axis_title=args.x_axis_title, x_axis_title_font_size=args.x_axis_title_font_size,
                     x_axis_tick_font_size=args.x_axis_tick_font_size, y_axis_tick_font_size=args.y_axis_tick_font_size)


if __name__ == "__main__":
    main()
