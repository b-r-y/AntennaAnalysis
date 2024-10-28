"""
Created on Nov 15, 2018

@author: bry
"""

import sys
import os
from typing import Any
import numpy as np
import unicodedata
import re
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.contour import QuadContourSet
from matplotlib.axes import Axes
from matplotlib.ticker import MaxNLocator
from .includes import utilities


class PatPlot:
    """
    Main plotting class.
    """

    def __init__(self, Patterns, legends):
        """
        This is the constructor for the plotting class.

        Arguments:
            fields (str): the type of fields pattern (list): list of ndarrays with fields data -\
                **MUST be in linear scale!!!** Use includes.patterns.prepare_for_plotting() to\
                generate the correct inputs.
            legends (list): a list of legend strings index matched to the ``patterns``

        Returns:
            object: returns an object of this class for plotting
        """
        self.patterns = Patterns
        self.legends = legends
        self.plot_settings = {
            "plot_type": "polar",
            "scale": "log",
            "dynamic_range": None,
            "plot_range": [[0, 180], [0, 360]],
            "output_path_pad": "",
        }

    @staticmethod
    def __apply_plot_range_limit(plot_range, pattern):
        """
        Warning:
            This is a private method and is NOT intended for external use.

        Applies a filter on the data so that all values outside the specified ``plot_range`` are \
        set to zero. This allows for plotting of only the area of interest.

        Arguments:
            plot_range ([[int,int],[int,int]]): specifies the plotting range of interest. All \
            values outside are set to zero. The values are in degrees. The two values pairs specify\
            the start and stop ranges for the Theta and Phi angles correspondingly.
            pattern (array): array with the fields data

        Returns:
            array: returns same size array with values outside the plot range set to zero
        """
        if plot_range is None:
            plot_range = [[0, 180], [0, 360]]

        th_step, ph_step = utilities.get_angle_steps(pattern, sphere_type="closed")

        selection_matrix = np.full(pattern.shape, np.finfo(float).eps)
        selection_matrix[
            int(plot_range[0][0] // th_step) : int(
                (plot_range[0][1] + th_step) // th_step
            ),
            int(plot_range[1][0] // ph_step) : int(
                (plot_range[1][1]) + ph_step // ph_step
            ),
        ] = 1
        pattern = selection_matrix * pattern

        return pattern

    @staticmethod
    def __process_data_select_inputs(data_select):
        """
        Converts the ``data_select`` inputs into info string for the title, ``plot_angle`` and
        ``plot_vector``.
        """
        if data_select.__class__ == str:
            if data_select == "XZ":
                info_str = data_select
                plot_angle = "theta"
                plot_vector = [0]
            elif data_select == "YZ":
                info_str = data_select
                plot_angle = "theta"
                plot_vector = [90]
            elif data_select == "XY":
                info_str = data_select
                plot_angle = "phi"
                plot_vector = [90]
            else:
                sys.exit(
                    "Undefined high level data_select option! --> {}".format(
                        data_select
                    )
                )
        else:
            plot_angle = data_select[0]
            plot_vector = data_select[1]
            if plot_angle == "theta":
                other_angle = "phi"
            elif plot_angle == "phi":
                other_angle = "theta"
            else:
                sys.exit("Unknown plot angle! --> {}".format(plot_angle))
            info_str = (
                "Sweeping "
                + plot_angle
                + " @ "
                + other_angle
                + " values"
                + plot_vector.__str__()
            )
        return info_str, plot_angle, plot_vector

    @staticmethod
    def __get_plot_data(pattern, plot_angle, curr_angle):
        """
        This method fetches the requested slices from the pattern to extract only the data needed
        for a specific plot.

        Arguments:
            pattern (array): the full 3D pattern with some field
            plot_angle (str): the angle to be swept for plotting - it can be theta or phi
            curr_angle (list): a list of angles for each of which the ``plot_angle`` will be swept

        Returns:
            array: returns a vector with the 360 degrees pattern slice. If the ``plot_angle`` is\
                theta it will wrap around the sphere and return the "other side" as well.
        """
        th_step, ph_step = utilities.get_angle_steps(pattern, sphere_type="closed")

        if plot_angle == "theta":
            ind_1 = int(round(curr_angle / ph_step))
            ind_2 = int(ind_1 + round(180 / ph_step))
            # Wrap around in case of angles for phi > 180 degrees
            if ind_2 > int(round(360 / ph_step)):
                ind_1 = ind_1 - int(round(360 / ph_step))
                ind_2 = ind_2 - int(round(360 / ph_step))
        elif plot_angle == "phi":
            ind_1 = int(round(curr_angle / th_step))
        else:
            sys.exit("Unknown ``plot_angle`` option!!! --> {}".format(plot_angle))

        if plot_angle == "theta":
            g_1 = pattern[:, ind_1]
            g_2 = pattern[:, ind_2]
            g_2 = np.flip(g_2, 0)
            pattern_slice = np.concatenate((g_1, g_2), axis=0)
        elif plot_angle == "phi":
            pattern_slice = pattern[ind_1, :].T  # /2 = 90 degrees

        return pattern_slice

    @staticmethod
    def __set_xaxis(axes, plot_angle, spacing=30, linear_plot_offset=0):
        """
        Set up the axis for standard 2D plots. Should not be used independently.
        """
        no_samples = int((180 / spacing) + 1)
        xticks = np.linspace(0, 360, int(2 * no_samples - 1)) * np.pi / 180
        if plot_angle == "theta":
            xlabels = [
                str(int(x))
                for x in np.linspace(0, 180, no_samples).tolist()
                + np.flip(
                    np.linspace(0, -180 + spacing, no_samples - 1), axis=0
                ).tolist()
            ]
            xlabels[180 // spacing] = r"$\pm 180$"
        elif plot_angle == "phi":
            xlabels = [
                str(int(x)) for x in np.linspace(0, 360, 2 * no_samples - 1).tolist()
            ]

        if axes.name == "polar":
            xticks = xticks[0:-1]
            xlabels = xlabels[0:-1]
        elif axes.name == "rectilinear":
            shift = int(linear_plot_offset / spacing)
            xlabels = np.roll(xlabels, shift).tolist()
            xlabels.insert(0, xlabels[-1])
            xlabels.pop(shift)

        axes.set_xticks(xticks)
        axes.set_xticklabels(xlabels)

    @staticmethod
    def __set_yaxis(axes, spacing=30, linear_plot_offset=0):
        """
        Set up teh Y axis mostly for the 2D contour plot.
        """
        no_ticks = int((180 / spacing) + 1)
        yticks = np.linspace(0, 180, int(no_ticks)) * np.pi / 180
        ylabels = [str(int(x)) for x in np.linspace(0, 180, no_ticks).tolist()]
        shift = int(linear_plot_offset / spacing)
        ylabels = np.roll(ylabels, shift).tolist()
        ylabels.insert(0, ylabels[-1])
        ylabels.pop(shift)

        axes.set_yticks(yticks)
        axes.set_yticklabels(ylabels)

    @staticmethod
    def __apply_linear_plot_offset(plot_axis, selected_data, linear_plot_offset):
        """
        Applies rolling offset in linear plots to bring areas of the curve into the center of the
        plot. Should not be used independently.
        """
        stitching_offset = plot_axis.__len__() - np.unique(plot_axis).__len__()
        resolution = 360 / (((plot_axis * 180 / np.pi).shape[0]) - 1 - stitching_offset)
        shift = int(linear_plot_offset / resolution + stitching_offset)
        selected_data = np.roll(selected_data, shift)
        return plot_axis, selected_data

    def __apply_dynamic_range_limit(self, plot_curves):
        """
        Warning:
            This is a private method and is NOT intended for external use.

        Applies the specified dynamic range limits and strategies. This methos does NOT modify any
        data but rather defines the plot axis limits.

        The method uses the object's internally stored dynamic range tuple. The tuple includes the\
        following definitions: [reference value strategy, [reference value / dynamic range, dynamic\
        range], scale]

            - Reference value strategy refers to the method by which the reference value for the\
                plot is determined (the maximum plot value). It can be one the following:

                    - 'auto' - when set to auto the input parameter ``plot_curves`` is used to scan\
                        the maximum value among the concrete context of curves. Requires a single\
                        value of range at the second place of the tuple.
                    - 'MaxFull3D' - this strategy will scan all fields in teh object and find the\
                        maximum value across the full pattern. Requires a single\
                        value of range at the second place of the tuple.
                    - 'manual' - value based specification for the reference. This requires two\
                        values list in the second place in teh tuple, where the first value is the\
                        manually defined reference and the second is the range.

            - the second place in the tuple can be either a single value defining the range in the\
                case of 'auto' or 'MaxFull3D' reference value strategy, or two values list defining\
                the reference value and range respectively, in the case of 'manual' reference value\
                strategy.
            -  the third value in the tuple can be either 'log' or 'lin' specifying logarithmic or\
                linear scale for the dynamic range definitions. Note that this is independent of \
                the actual scale for the plot.

        Arguments:
            plot_curves (tuple): a list of multiple plot curves. This parameter is used when the\
                ``auto`` option for dynamic range determination is used.

        Returns:
            list: returns a list with two numeric value with the MIN and MAX for teh plot axis \
                range

        Examples:
            dynamic_range = ['auto', 40, 'log'] is the default if nothing else is specified. This \
            is interpreted as 40 dB dynamic range automatically determined based on the curves \
            being plotted. Note that if a specific radaition pattern slice of particularly low \
            power is being plotted this plot maximum can be 10-20 dB lower than the main gain\
            direction not captured in the slice.

            dynamic_range = ['auto', 10000, 'lin'] is equivalent to the default in linear scale.

            dynamic_range = ['manual', [10, 30], 'log'] defines manual maximum for the plot of 10 dB
            and minimum for the plot 30 dB below at -20 dB. Note that in this definitions a peak
            gain above 10 dB might be out of the plot.

            In any of the above examples the actual plot can be in linear or logarithmic. The
            'lin'/'log' specified here only refers to the scale for the dynamic range definition.

        """
        if self.plot_settings.get("dynamic_range") is None:
            self.plot_settings["dynamic_range"] = ["auto", 40, "log"]

        # First handle the scale
        # expects a single value for the dynamic range
        if self.plot_settings.get("dynamic_range")[0] == "auto":
            lower_limit = self.plot_settings.get("dynamic_range")[1]
            # following line can be log or lin depending on input
            ref_level = np.max([curr_p.max() for curr_p in plot_curves])
        # expects a reference level and dynamic range
        elif self.plot_settings.get("dynamic_range")[0] == "MaxFull3D":
            lower_limit = self.plot_settings.get("dynamic_range")[1]
            ref_level = np.max([curr_p.max() for curr_p in self.patterns])
        # expects a reference level and dynamic range
        elif self.plot_settings.get("dynamic_range")[0] == "manual":
            lower_limit = self.plot_settings.get("dynamic_range")[1][1]
            ref_level = float(self.plot_settings.get("dynamic_range")[1][0])
        else:
            sys.exit("Unknown dynamic range option")

        if self.plot_settings.get("dynamic_range")[2] == "log":
            lower_limit = 10 ** (lower_limit / 10)
            if self.plot_settings.get("dynamic_range")[0] == "manual":
                ref_level = 10 ** (ref_level / 10)
        elif self.plot_settings.get("dynamic_range")[2] == "lin":
            pass
        else:
            sys.exit("Unknown dynamic range option")

        if not self.plot_settings.get("dynamic_range")[0] == "manual":
            plot_max = ref_level * 1.1
        else:
            plot_max = ref_level
        plot_min = ref_level / lower_limit

        plot_limits = [plot_min, plot_max]

        if self.plot_settings.get("scale") == "log":
            plot_limits = 10 * np.log10(plot_limits)

        plot_limits[0] = np.floor(plot_limits[0])
        plot_limits[1] = np.ceil(plot_limits[1])

        return plot_limits

    def __apply_scale(self, field_index, field_data):
        """
        Warning:
            This is a private method and is NOT intended for external use.

        Converts the data to be plotted to logarithmic scale if this is specified in the object's
        ``scale`` parameter. If the requested plot is to be in linear scale this function does
        nothing as the data is already assumed to be supplied in linear scale.

        Arguments:
            field_data (array): array with the pattern data.

        Returns:
            array: returns the same array converted to logarithmic scale.
        """
        if self.plot_settings.get("scale") == "log":
            field_data = utilities.lin2log(self.legends[field_index][0], field_data)
        return field_data

    def __make_title(self, axes, info_str):
        """
        Creates the figure title.
        """
        if self.plot_settings.get("scale") == "log":
            scale = "dB"
        else:
            scale = self.plot_settings.get("scale")
        title_pad = 10
        info_str = info_str + " [" + scale + "]"
        title = axes.set_title(info_str, pad=title_pad)
        return info_str, title

    def __make_axes(self, fig, plot_angle="theta", override=None):
        """
        Makes the necessary axes for plotting. In case of polar axes the view is always from the
        positive direction of the normal axis according to the IEEE definitions.
        """
        if override is None:
            if self.plot_settings.get("plot_type") == "polar":
                projection = self.plot_settings.get("plot_type")
            else:
                projection = "rectilinear"
            ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.9], projection=projection)

            if self.plot_settings.get("plot_type") == "polar":
                if plot_angle == "phi":
                    ax1.set_theta_zero_location("S")
                elif plot_angle == "theta":
                    ax1.set_theta_zero_location("N")
        elif override == "stats":
            ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.9], projection="rectilinear")
        elif override == "contour":
            ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.9], projection="rectilinear")
        elif override == "3d":
            ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.9], projection="3d")
            # ax1 = Axes3D(fig,rect = [0.1,0.1,0.8,0.9])
        return ax1

    def __get_standard_plot_axis(self, plot_angle):
        """
        Generates the plot axis for a 2D plot depending on the angle selected.
        """
        # Extracts the shape of the plotting object's patterns. Since the object is always
        # initialized with some patterns loaded and those MUST be of equal size this function simply
        # gets that from the first one.
        no_samples = self.patterns[0].shape

        if plot_angle == "theta":
            plot_axis = (
                np.concatenate(
                    (
                        np.linspace(0, 180, no_samples[0]),
                        np.linspace(0, 180, no_samples[0]) + 180,
                    ),
                    axis=0,
                )
                * np.pi
                / 180
            )
        elif plot_angle == "phi":
            plot_axis = np.linspace(0, 360, no_samples[1]) * np.pi / 180
        else:
            sys.exit("Unknown ``plot_angle`` option!!! --> {}".format(plot_angle))

        return plot_axis

    def __augment_legends(self, plot_angle, plot_vector):
        """
        Augments legends with specifics of teh plot slices used.
        """
        plot_legends = []

        if plot_angle == "theta":
            other_anlge = "phi"
        else:
            other_anlge = "theta"

        for i_0 in range(self.legends.__len__()):
            for i_1 in range(plot_vector.__len__()):
                plot_legends.append(
                    "["
                    + other_anlge
                    + "="
                    + str(plot_vector[i_1])
                    + "deg] --> "
                    + " ".join(self.legends[i_0])
                )

        return plot_legends

    def __plot_standard_2d(self, linear_plot_offset, plot_angle, plot_vector, ax1):
        """
        Plots the standard polar or linear plots.
        """
        data_to_plot = []
        plot_axis = self.__get_standard_plot_axis(plot_angle)
        for idx, curr_p in enumerate(self.patterns):
            curr_p = self.__apply_plot_range_limit(
                self.plot_settings.get("plot_range"), curr_p
            )
            for curr_angle in plot_vector:
                selected_data = self.__get_plot_data(curr_p, plot_angle, curr_angle)
                if self.plot_settings.get("plot_type") == "linear":
                    plot_axis, selected_data = self.__apply_linear_plot_offset(
                        plot_axis, selected_data, linear_plot_offset
                    )
                data_to_plot.append(selected_data)
                # the scale application is BEFORE the plotting so that the correct scale is plotted.
                # Yet it is AFTER the storing into the data list so that the dynamic range
                # is applied on linear data!
                selected_data = self.__apply_scale(idx, selected_data)
                ax1.plot(plot_axis, selected_data)

        y_plot_limits = self.__apply_dynamic_range_limit(data_to_plot)

        return y_plot_limits

    def __plot_2d_contour(self, ax1, ind, n_bins, linear_plot_offset):
        """
        Plots 2D contour via ``contourf``
        """
        curr_p = self.patterns[ind]
        curr_p = self.__apply_plot_range_limit(
            self.plot_settings.get("plot_range"), curr_p
        )
        data_to_plot = self.__apply_scale(ind, curr_p)
        data_to_plot = np.roll(data_to_plot, linear_plot_offset, axis=(0, 1))
        z_limit = self.__apply_dynamic_range_limit(curr_p)
        cmap = plt.get_cmap("jet")
        levels = MaxNLocator(nbins=n_bins).tick_values(z_limit[0], z_limit[1])
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        vth, vph = utilities.get_angle_vectors(data_to_plot, "closed")
        mesh = ax1.contourf(vph.T, vth.T, data_to_plot, levels, cmap=cmap, norm=norm)
        return mesh

    @staticmethod
    def __spherical_to_cartesian(data):
        """
        Converts from spherical to cartesian coordinates.
        """
        vth, vph = utilities.get_angle_vectors(data, "closed")

        x_data = (data / data.max()) * np.sin(vth.T) * np.cos(vph.T)
        y_data = (data / data.max()) * np.sin(vth.T) * np.sin(vph.T)
        z_data = (data / data.max()) * np.cos(vth.T)

        return x_data, y_data, z_data

    def __plot_3d(self, ax1, data_to_plot, cmap):
        """
        Converts from spherical to cartesian coordinates and plots a ``plot_surface`` with patch
        faces painted according to ``cmap``.
        """
        x_data, y_data, z_data = self.__spherical_to_cartesian(data_to_plot)
        ax1.plot_surface(
            x_data, y_data, z_data, facecolors=cmap(data_to_plot / data_to_plot.max())
        )

    def __save_figure(self, figure_file_name, fig, extra_artists: list = None):
        """
        Saves the figures in PNG, PDF and SVG file formats.
        """
        save_folder = "output/" + self.plot_settings.get("output_path_pad")
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        figure_file_name = self.__slugify(figure_file_name, True)
        fig.savefig(
            save_folder + "/" + figure_file_name + ".png",
            bbox_extra_artists=extra_artists,
            bbox_inches="tight",
        )
        fig.savefig(
            save_folder + "/" + figure_file_name + ".pdf",
            bbox_extra_artists=extra_artists,
            bbox_inches="tight",
        )
        fig.savefig(
            save_folder + "/" + figure_file_name + ".svg",
            bbox_extra_artists=extra_artists,
            bbox_inches="tight",
        )

    def set_plot_params(self, plot_settings: dict[str, Any]):
        """
        Sets the plot parameters for this plot object. All plots of the same object share these
        parameters. It is possible to create a single PatPlot object and alternate plots and
        settings updates to create different plots from the same data.

        Arguments:
            plot_settings (dict): A dictionary with plot settings keys and values. Only the ones\
                specified will be changed - the rest remain unaltered. The possible keys/values\
                are:

                - **'plot_type'** - *'polar'* (Default) and 'linear'. Defines polar or linear plots\
                        for the standard 2D plots. This is ignored for 2D or 3D plots.

                - **'scale'** - *'log'* (Default) and 'lin'. Defines logarithmic or linear data \
                    plot along the Y axis.

                - **'dynamic_range'** - *None* (Default). See ``__apply_dynamic_range_limit`` for\
                    details on how to define the dynamic range.

                - **'plot_range'** - *[[0,180],[0,360]]* (Default).See ``__apply_plot_range_limit``\
                     for details on how to define the plotting range.

                - **'output_path_pad'** - ''(Default). A string defining additional subfoldering in\
                    the output folder.

        Returns:
            int: updates the internal settings to be used and returns 0 if successful.
        """
        self.plot_settings = utilities.update_params(self.plot_settings, plot_settings)

        return 0

    @staticmethod
    def __get_wraparound_axis_limits(pos, neg):
        pos[pos > 360] = pos[pos > 360] - 360
        neg[neg >= 360] = neg[neg >= 360] - 360

        lim = np.array([pos, neg])
        lim = np.array([np.min(lim), np.max(lim)])
        lim.sort()

        return lim

    @staticmethod
    def __slugify(value, allow_unicode=False):
        """
        Taken from https://github.com/django/django/blob/master/django/utils/text.py
        Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
        dashes to single dashes. Remove characters that aren't alphanumerics,
        underscores, or hyphens. Convert to lowercase. Also strip leading and
        trailing whitespace, dashes, and underscores.

        To account for the frequency as a float dots are also allowed.
        """
        value = str(value)
        if allow_unicode:
            value = unicodedata.normalize("NFKC", value)
        else:
            value = (
                unicodedata.normalize("NFKD", value)
                .encode("ascii", "ignore")
                .decode("ascii")
            )
        value = re.sub(r"[^\w\.\s-]", "", value.lower())
        return re.sub(r"[-\s]+", "-", value).strip("-_")

    def __set_axis_limits(self, plot_angle, linear_plot_offset):
        if self.plot_settings.get("plot_type") == "linear":
            if plot_angle == "theta":
                lim = self.__get_wraparound_axis_limits(
                    np.array(self.plot_settings.get("plot_range")[0])
                    + linear_plot_offset,
                    360
                    - np.array(self.plot_settings.get("plot_range")[0])
                    + linear_plot_offset,
                )
            elif plot_angle == "phi":
                lim = (
                    np.array(self.plot_settings.get("plot_range")[1])
                    + linear_plot_offset
                )
                if all(lim >= 360):
                    lim = lim - 360
                elif any(lim > 360):
                    pos = np.array([np.min(lim), 360])
                    neg = np.array([0, np.max(lim) - 360])
                    lim = self.__get_wraparound_axis_limits(pos, neg)
            lim = lim * np.pi / 180
        elif self.plot_settings.get("plot_type") == "polar":
            lim = np.array([0, 360]) * np.pi / 180
        return lim

    def make_standard_plot(
        self,
        data_select="XZ",
        xtick_spacing=30,
        linear_plot_offset=0,
        customizations: dict = None,
    ) -> object:
        """
        Standard plotting method for the PatPLot object - makes plots mostly according to the IEEE
        standard definitions for antennas and primary planes.

        Arguments:
            data_select (tuple): describes the slice of the full sphere to be plotted. The first\
                input is the sweep angle - possible values are theta and phi. The second input is\
                the non-sweep/'other' angle values in degrees at which values the sweep is\
                performed. In the case of theta sweep angle both +theta and -theta sides are taken.

                Short hand notations for the IEEE primary planes can also be defined as strings -\
                XZ, YZ and XY. These are equivalent to the following full definitions:

                - 'XZ' = ['theta', [0]] - sweep angle theta at phi 0 degrees.
                - 'YZ' = ['theta', [90]] - sweep angle theta at phi 90 degrees.
                - 'XY' = ['phi', [90]] - sweep angle phi at theta 90 degrees.

            xtick_spacing (int): defines the x axis major ticks. Since these plots are angles this\
                effectively sets the angular resolution of the grid of the plot.

            linear_plot_offset (int): defines the degrees of offset to be applied to the data along\
                the **X axis** in case of linear plots. This is particular useful for linear plots\
                of gain in case the antenna primary gain is in the theta=0 degrees direction -\
                typical patch antenna. Setting this to 180 degrees brings the peak at the center \
                of the plot.

            customizations (dict): a dictionary specifying various plot customizations typically\
                used when generating plots for datasheets or customer queries. See \
                :meth:`~includes.plotting.PatPlot._parse_customizations` for availabe options.
        Returns:
            object: returns a handle to the newly created figure. It also outputs the figure in\
                PDF/SVG and PNG formats at the output folder in the settings.
        """
        info_str, plot_angle, plot_vector = self.__process_data_select_inputs(
            data_select
        )
        fig = plt.figure()
        ax1 = self.__make_axes(fig, plot_angle)
        info_str, _ = self.__make_title(ax1, info_str)

        y_plot_limits = self.__plot_standard_2d(
            linear_plot_offset, plot_angle, plot_vector, ax1
        )

        self.__set_xaxis(ax1, plot_angle, xtick_spacing, linear_plot_offset)

        ax1.margins(x=0, y=0)
        ax1.set_xlim(self.__set_axis_limits(plot_angle, linear_plot_offset))
        ax1.set_ylim(y_plot_limits)

        lgd = self.__augment_legends(plot_angle, plot_vector)
        lgd = ax1.legend(
            lgd,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.05),
            fancybox=True,
            shadow=True,
            ncol=1,
        )
        box = ax1.get_position()
        ax1.set_position(
            [box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.8]
        )

        plt.grid(True)
        extra_artists = [lgd]

        if customizations is not None:
            self._parse_customizations(ax1, customizations)
            if (
                "show legend" in customizations.keys()
                and customizations["show legend"] is False
            ):
                extra_artists = None

        self.__save_figure(
            info_str + " " + self.plot_settings.get("plot_type"), fig, extra_artists
        )

        plt.show()

        return fig

    def _parse_customizations(self, axes: object, customizations: dict = None) -> None:
        """
        Allows a number of plot customizations to be passed as a dictionary and augment \
        the default plotting behavior.

        Arguments:
            axes (object): the axes object of the plot. In the case of contour plots the\
                `QuadContourSet` return class should be passed to allow color map and \
                normalization changes.

            customizations (dict): a dictionary defining the parameters to modify and \
                its new value. Availabe options are:

                - {'title': 'str'} - sets the title to `str`
                - {'xlabel': 'str'} - sets the primary X axis label to `str`
                - {'ylabel': 'str'} - sets the primary Y axis label to `str`
                - {'show legend': bool} - by default the legend is shown with a lot of\
                    engineering detail. Setting this false will hide it.
                - {'set Colormap': object} - an object of class `matplotlib.colors.Colormap\
                    <https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.Colormap.html#matplotlib.colors.Colormap>`_ \
                    to be used.
                - {'set Normalization': object} - an object of one of the normalization classes in `matplotlib.colors\
                    <https://matplotlib.org/stable/api/colors_api.html>`_ to be used.

        """
        if type(customizations) == dict:
            for k in customizations.keys():
                if isinstance(axes, Axes):
                    if k == "ylabel":
                        if self.plot_settings.get("plot_type") == "polar":
                            label_offset = 20
                        else:
                            label_offset = None
                        axes.set_ylabel(customizations[k], labelpad=label_offset)
                    elif k == "xlabel":
                        axes.set_xlabel(customizations[k])
                    elif k == "title":
                        axes.set_title(customizations[k])
                    elif k == "show legend":
                        if not customizations[k]:
                            axes.get_legend().remove()
                if isinstance(axes, QuadContourSet):
                    if k == "set Colormap":
                        axes.set_cmap(customizations[k])
                    elif k == "set Normalization":
                        axes.set_norm(customizations[k])

    def make_statistics_plot(self, stat_plot_type="CDF", customizations: dict = None):
        """Make a CDF or CCDF plots fo the data in the PatPlot object.

        Arguments:
            stat_plot_type (str): type of the plot. Possible values are \
                'CDF' (default) and 'CCDF'

            customizations (dict): a dictionary specifying various plot \
                customizations typically used when generating plots for \
                datasheets or customer queries. See \
                :meth:`~antenna_analysis.plotting.PatPlot._parse_customizations` \
                for availabe options.

        Returns:
            object: returns a handle to the newly created figure. It also outputs the figure in\
                PDF/SVG and PNG formats at the output folder in the settings.
        """
        fig = plt.figure()
        ax1 = self.__make_axes(fig, override="stats")
        info_str, _ = self.__make_title(ax1, stat_plot_type)

        data_to_plot = []
        for idx, curr_p in enumerate(self.patterns):
            curr_p = utilities.open_sphere(curr_p)
            curr_p = utilities.apply_analysis_range(
                curr_p, self.plot_settings.get("plot_range")
            )
            sorted_data = np.sort(curr_p)
            data_to_plot.append(sorted_data)
            sorted_data = self.__apply_scale(idx, sorted_data)
            if stat_plot_type == "CDF":
                yvals = np.arange(len(sorted_data)) / float(len(sorted_data) - 1)
            elif stat_plot_type == "CCDF":
                yvals = 1 - (np.arange(len(sorted_data)) / float(len(sorted_data) - 1))
            ax1.plot(sorted_data, yvals)

        ax1.set_xlim(self.__apply_dynamic_range_limit(data_to_plot))
        ax1.set_ylim((0, 1))

        lgd = []
        for i_0 in range(self.legends.__len__()):
            lgd.append(" ".join(self.legends[i_0]))

        lgd = ax1.legend(
            lgd,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.05),
            fancybox=True,
            shadow=True,
            ncol=1,
        )
        box = ax1.get_position()
        ax1.set_position(
            [box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.8]
        )

        plt.minorticks_on()
        plt.grid(visible=True, which="both")

        if customizations is not None:
            self._parse_customizations(ax1, customizations)
            if (
                "show legend" in customizations.keys()
                and customizations["show legend"] is False
            ):
                extra_artists = None

        plt.show(block=True)

        self.__save_figure(info_str, fig, [lgd])

        return fig

    def make_countour_plot(
        self,
        ind,
        n_bins=30,
        tick_spacing=None,
        linear_plot_offset=None,
        customizations: dict = None,
    ) -> object:
        """
        Makes a 2D contour plot from the plotting object at index ``ind`` and using ``n_bins``
        color map levels.

        Note:
            the single index used for this function refers to the fields the plotting object is\
            created with and NOT the full electric fields index from the radiation patterns object.

        Arguments:
            ind (int): index of the field for analysis of the plotting object

            n_bins (int): number of bins for the colormap. Default is 30.

            tick_spacing ([int,int]): defines the theta and phi tick spacing. Default is 30 degrees
                for both angles.

            linear_plot_offset ([int,int]): defines the wrapping around of the plot along the theta\
                and phi angles. It can be used to center a main radiation direction. Default values\
                are zero degrees shift for both angles.

            customizations (dict): a dictionary specifying various plot customizations typically\
                used when generating plots for datasheets or customer queries. See \
                :meth:`~antenna_analysis.plotting.PatPlot._parse_customizations`\
                for availabe options.

        Returns:
            object: returns a handle to the newly created figure. It also outputs the figure in\
                PDF/SVG and PNG formats at the output folder in the settings.
        """
        if tick_spacing is None:
            tick_spacing = [30, 30]
        if linear_plot_offset is None:
            linear_plot_offset = [0, 0]

        fig = plt.figure()
        ax1 = self.__make_axes(fig, override="contour")
        info_str, _ = self.__make_title(ax1, "2D Contour")

        mesh = self.__plot_2d_contour(ax1, ind, n_bins, linear_plot_offset)

        plt.xlabel("Phi [deg]", labelpad=0)
        plt.ylabel("Theta [deg]", labelpad=0)
        plt.gca().invert_yaxis()
        self.__set_xaxis(ax1, "phi", tick_spacing[1], linear_plot_offset[1])
        self.__set_yaxis(ax1, tick_spacing[0], linear_plot_offset[0])

        lgd = ax1.legend(
            [plt.Rectangle((0, 0), 1, 1)],
            [" ".join(self.legends[ind])],
            loc="upper center",
            bbox_to_anchor=(0.5, -0.10),
            fancybox=True,
            shadow=True,
            ncol=1,
        )
        box = ax1.get_position()
        ax1.set_position(
            [box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.8]
        )
        fig.colorbar(mesh, ax=ax1)
        extra_artists = [lgd]
        plt.grid(visible=False, which="both")

        if customizations is not None:
            self._parse_customizations(mesh, customizations)
            self._parse_customizations(ax1, customizations)
            if (
                "show legend" in customizations.keys()
                and customizations["show legend"] is False
            ):
                extra_artists = None

        plt.show(block=True)

        self.__save_figure(
            info_str + " -> " + "".join(self.legends[ind]), fig, extra_artists
        )

        return fig

    def make_3d_plot(self, ind, n_bins=30):
        """
        Makes a 3D pattern plot from the plotting object at index ``ind`` and using ``n_bins``
        color map levels.

        Note:
            the single index used for this function refers to the fields the plotting object is\
            created with and NOT the full electric fields index from the radiation patterns object.

        Arguments:
            ind (int): index of the field for analysis of the plotting object
            n_bins (int): number of bins for the colormap. Default is 30.

        Returns:
            object: returns a handle to the newly created figure. It also outputs the figure in\
                PDF/SVG and PNG formats at the output folder in the settings.
        """
        fig = plt.figure()
        ax1 = self.__make_axes(fig, override="3d")
        info_str, _ = self.__make_title(ax1, "3D Sphere")

        curr_p = self.patterns[ind]
        curr_p = self.__apply_plot_range_limit(
            self.plot_settings.get("plot_range"), curr_p
        )
        data_to_plot = self.__apply_scale(ind, curr_p)
        z_limit = self.__apply_dynamic_range_limit(curr_p)
        data_to_plot[data_to_plot < z_limit[0]] = z_limit[0]
        data_to_plot = data_to_plot + np.abs(z_limit[0])

        cmap = plt.get_cmap("jet", n_bins)
        color_map = cm.ScalarMappable(cmap=cmap)
        color_map.set_array(data_to_plot - np.abs(z_limit[0]))

        self.__plot_3d(ax1, data_to_plot, cmap)

        # ax1.set_aspect('equal') # removed as not implemented yet. Older version produced wrong results... see https://github.com/matplotlib/matplotlib/issues/15382
        ax1.axes.set_xlim([-1, 1])
        ax1.axes.set_ylim([-1, 1])
        ax1.axes.set_zlim([-1, 1])
        ax1.axis("off")

        lgd = ax1.legend(
            [plt.Rectangle((0, 0), 1, 1)],
            [" ".join(self.legends[ind])],
            loc="upper center",
            bbox_to_anchor=(0.5, -0.10),
            fancybox=True,
            shadow=True,
            ncol=1,
        )
        box = ax1.get_position()
        ax1.set_position(
            [box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.8]
        )
        fig.colorbar(color_map, ax=ax1)

        plt.grid(visible=False, which="both")
        plt.show(block=True)

        self.__save_figure(info_str, fig, [lgd])

        return fig
