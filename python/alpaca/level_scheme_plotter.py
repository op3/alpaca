# This file is part of alpaca.

# alpaca is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# alpaca is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with alpaca.  If not, see <https://www.gnu.org/licenses/>.

# Copyright (C) 2021 Udo Friman-Gayer

import matplotlib.pyplot as plt
import numpy as np


class LevelSchemePlotter:
    """Class to plot a labeled level scheme with a single excitation and a decay cascade

    This class is intended to visualize the transition cascades whose angular correlations are
    given by the AngularCorrelation class.
    The first transition, which defines the 'orientation' of the system, is assumed to be an
    excitation from the ground state.
    All other transitions belong to a cascade whose last transition is assumed to be observed.

    The visualization is realized using the matplotlib package.
    LevelSchemePlotter draws the level scheme on a matplotlib.axes.Axes object that is passed to its __init__ function.

    A plotted level scheme contains the following elements:

    - Horizontal lines that represent states.
    - State labels that indicate the spin- and parity quantum numbers.
    - One upward arrow that represents the excitation.
    - Downward arrows that represent the decays.
    - Transition labels that indicate the EM character and the multipolarity.
    - Transition labels that indicate the Multipole-mixing ratio of that transition.

    Attributes
    ----------

    ax: matplotlib.axes.Axes object
        Axes on which the level scheme should be drawn
    ini_sta: State
        Initial state of the cascade. Same format as the corresponding input for the AngularCorrelation class.
    cas_ste: array of [Transition, State] pairs
        Cascade steps, given as a list of arbitrary length which contains Transition-State pairs.
        The first and the last transition of this list are assumed to be observed.
        Same format as the corresponding input for the AngularCorrelation class.
    min_x, max_x, min_y, max_y: float
        Limits of the x - and y axis as given by ax.
    range_x, range_y:
        Range of the x - and y axis.
    del_lab: list of str
        Labels for the multipole mixing ratios. The length of the list should equal the number of cascade steps.
    returns_to_initial_state: bool
        Determines whether the last step of the cascade should go back to the ground state.
        This option allows to make the distinction between a ground-state decay and a cascade that
        ends up in a state with the same quantum numbers as the ground state. (default: False)
    show_polarization: list of bool
        Determines which transition labels should indicate a polarization.
        This allows to indicate for which transition polarization information is available.
    fontsize: int or float
        Font size of the figure. Several other font sizes are scaled to this value.
    fontsize_single_multipole: float
        Font size for transition labels with a single multipole (default: 0.9*fontsize).
    fontsize_single_multipole: float
        Font size for transition labels with two possible multipoles. (default: 1.2*fontsize)
    """

    def __init__(
        self,
        axis,
        initial_state,
        cascade_steps,
        delta_labels,
        returns_to_initial_state=False,
        show_polarization=None,
        fontsize=12,
        state_line_width=2,
        arrow_width=2,
        offset=(0, 0),
        transition_label_rotation=90,
        em_variable_symbol="\sigma",
        parity_variable_symbol="\pm",
    ):
        self.ax = axis
        self.min_x, self.max_x = axis.get_xlim()
        self.range_x = self.max_x - self.min_x
        self.min_y, self.max_y = axis.get_ylim()
        self.range_y = self.max_y - self.min_y
        self.ini_sta = initial_state
        self.cas_ste = cascade_steps
        self.del_lab = delta_labels
        self.returns_to_initial_state = returns_to_initial_state
        self.show_polarization = [False] * len(cascade_steps)
        if show_polarization is not None:
            self.show_polarization = show_polarization

        ## Parameters for the plot
        # Fonts
        self.fontsize = fontsize
        self.fontsize_single_multipole = 0.9 * fontsize
        self.fontsize_two_multipoles = 1.2 * fontsize

        # State lines
        self.state_line_width = state_line_width
        self.state_x = 0.4 * self.range_x + self.min_x + offset[0] * self.range_x
        self.state_width = 0.4 * self.range_x
        self.intermediate_state_x = (
            0.55 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.intermediate_state_width = 0.25 * self.range_x

        self.initial_state_y = (
            0.2 * self.range_y + self.min_y + offset[1] * self.range_y
        )
        self.excited_state_y = (
            0.8 * self.range_y + self.min_y + offset[1] * self.range_y
        )

        # State labels
        self.state_label_left_x = (
            0.25 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.state_label_right_x = (
            0.9 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.parity_variable_symbol = parity_variable_symbol

        # Transition arrows
        self.arrow_width = arrow_width
        self.excitation_arrow_x = (
            0.5 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.decay_arrow_x = 0.7 * self.range_x + self.min_x + offset[0] * self.range_x
        self.arrow_head_length = 0.04 * self.range_y
        self.arrow_head_width = 0.03 * self.arrow_width
        self.excitation_arrow_color = "blue"
        self.decay_arrow_color = "red"

        # Transition labels
        self.em_variable_symbol = em_variable_symbol
        self.decay_label_right_x = (
            0.85 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.excitation_label_left_x = (
            0.18 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.transition_label_rotation = transition_label_rotation

        # Multipole mixing ratio (delta) labels
        self.delta_label_left_x = (
            0.4 * self.range_x + self.min_x + offset[0] * self.range_x
        )
        self.delta_label_right_x = (
            0.74 * self.range_x + self.min_x + offset[0] * self.range_x
        )

        # Order of drawing
        self.zorder_states = 0
        self.zorder_arrows = 1

    def plot(self):
        # Initial and excited state
        self.ax.plot(
            [self.state_x, self.state_x + self.state_width],
            [self.initial_state_y] * 2,
            color="black",
            linewidth=self.state_line_width,
            zorder=self.zorder_states,
        )
        self.ax.text(
            self.state_label_left_x,
            self.initial_state_y,
            self.ini_sta.tex(parity_variable_symbol=self.parity_variable_symbol),
            verticalalignment="center",
            fontsize=self.fontsize,
        )
        self.ax.plot(
            [self.state_x, self.state_x + self.state_width],
            [self.excited_state_y] * 2,
            color="black",
            linewidth=self.state_line_width,
            zorder=self.zorder_states,
        )
        self.ax.text(
            self.state_label_left_x,
            self.excited_state_y,
            self.cas_ste[0][1].tex(parity_variable_symbol=self.parity_variable_symbol),
            verticalalignment="center",
            fontsize=self.fontsize,
        )

        # Excitation
        self.ax.plot(
            [self.excitation_arrow_x] * 2,
            [self.initial_state_y, self.excited_state_y - self.arrow_head_length],
            "-",
            linewidth=self.arrow_width,
            color=self.excitation_arrow_color,
            zorder=self.zorder_arrows,
        )
        self.ax.arrow(
            self.excitation_arrow_x,
            self.excited_state_y - self.arrow_head_length,
            0.0,
            1.0e-5 * self.range_y,
            head_length=self.arrow_head_length,
            head_width=self.arrow_head_width,
            facecolor=self.excitation_arrow_color,
            edgecolor=self.excitation_arrow_color,
            zorder=self.zorder_arrows,
        )
        self.ax.text(
            self.excitation_label_left_x,
            0.5 * (self.excited_state_y - self.initial_state_y) + self.initial_state_y,
            self.cas_ste[0][0].tex(
                em_variable_symbol=self.em_variable_symbol,
                always_show_secondary=False,
                show_polarization=self.show_polarization[0],
            ),
            verticalalignment="center",
            fontsize=self.fontsize_single_multipole
            if self.cas_ste[0][0].delta == 0.0
            else self.fontsize_two_multipoles,
        )
        self.ax.text(
            self.delta_label_left_x,
            0.5 * (self.excited_state_y - self.initial_state_y) + self.initial_state_y,
            self.del_lab[0],
            verticalalignment="center",
            fontsize=self.fontsize,
            rotation=self.transition_label_rotation,
        )

        # Calculate position of states in decay cascade
        n_decay_steps = len(self.cas_ste) - 1
        excitation_delta_y = self.excited_state_y - self.initial_state_y
        cascade_states_y = np.arange(1, n_decay_steps + 1)

        if self.returns_to_initial_state:
            cascade_states_y = (
                self.excited_state_y
                - cascade_states_y * excitation_delta_y / (n_decay_steps)
            )
        else:
            cascade_states_y = (
                self.excited_state_y
                - cascade_states_y * excitation_delta_y / (n_decay_steps + 1)
            )

        # States in cascade
        for i in range(
            n_decay_steps if not self.returns_to_initial_state else n_decay_steps - 1
        ):
            self.ax.plot(
                [
                    self.intermediate_state_x,
                    self.intermediate_state_x + self.intermediate_state_width,
                ],
                [cascade_states_y[i]] * 2,
                color="black",
                linewidth=self.state_line_width,
                zorder=self.zorder_states,
            )
            self.ax.text(
                self.state_label_right_x,
                cascade_states_y[i],
                self.cas_ste[i + 1][1].tex(
                    parity_variable_symbol=self.parity_variable_symbol
                ),
                verticalalignment="center",
                fontsize=self.fontsize,
            )

        # First transition in cascade
        self.ax.plot(
            [self.decay_arrow_x] * 2,
            [self.excited_state_y, cascade_states_y[0] + self.arrow_head_length],
            "--" if n_decay_steps > 1 else "-",
            linewidth=self.arrow_width,
            color=self.decay_arrow_color,
            zorder=self.zorder_arrows,
        )
        self.ax.arrow(
            self.decay_arrow_x,
            cascade_states_y[0] + self.arrow_head_length,
            0.0,
            -1.0e-5 * self.range_y,
            head_length=self.arrow_head_length,
            head_width=self.arrow_head_width,
            color=self.decay_arrow_color,
            zorder=self.zorder_arrows,
        )
        self.ax.text(
            self.decay_label_right_x,
            0.5 * (self.excited_state_y - cascade_states_y[0]) + cascade_states_y[0],
            self.cas_ste[1][0].tex(
                em_variable_symbol=self.em_variable_symbol,
                always_show_secondary=False,
                show_polarization=self.show_polarization[1],
            ),
            verticalalignment="center",
            fontsize=self.fontsize_single_multipole
            if self.cas_ste[1][0].delta == 0.0
            else self.fontsize_two_multipoles,
        )
        self.ax.text(
            self.delta_label_right_x,
            0.5 * (self.excited_state_y - cascade_states_y[0]) + cascade_states_y[0],
            self.del_lab[1],
            verticalalignment="center",
            fontsize=self.fontsize,
            rotation=self.transition_label_rotation,
        )

        # Transitions in cascade
        for i in range(1, n_decay_steps):
            self.ax.plot(
                [self.decay_arrow_x] * 2,
                [cascade_states_y[i - 1], cascade_states_y[i] + self.arrow_head_length],
                "--" if i < n_decay_steps - 1 else "-",
                linewidth=self.arrow_width,
                color=self.decay_arrow_color,
                zorder=self.zorder_arrows,
            )
            self.ax.arrow(
                self.decay_arrow_x,
                cascade_states_y[i] + self.arrow_head_length,
                0.0,
                -1.0e-5 * self.range_y,
                head_length=self.arrow_head_length,
                head_width=self.arrow_head_width,
                color=self.decay_arrow_color,
                zorder=self.zorder_arrows,
            )
            self.ax.text(
                self.decay_label_right_x,
                0.5 * (cascade_states_y[i - 1] - cascade_states_y[i])
                + cascade_states_y[i],
                self.cas_ste[i + 1][0].tex(
                    em_variable_symbol=self.em_variable_symbol,
                    always_show_secondary=False,
                    show_polarization=self.show_polarization[i + 1],
                ),
                verticalalignment="center",
                fontsize=self.fontsize_single_multipole
                if self.cas_ste[i + 1][0].delta == 0.0
                else self.fontsize_two_multipoles,
            )
            self.ax.text(
                self.delta_label_right_x,
                0.5 * (self.excited_state_y - cascade_states_y[i - 1])
                + cascade_states_y[i],
                self.del_lab[i + 1],
                verticalalignment="center",
                fontsize=self.fontsize,
                rotation=self.transition_label_rotation,
            )
