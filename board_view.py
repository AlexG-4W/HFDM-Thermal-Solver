"""
Board canvas: owns the Matplotlib figure and every artist drawn on it.

Why this module exists
----------------------
Rendering used to be a side effect of user actions. Four methods in gui.py -
view_topology(), on_table_edit(), start_simulation() and update_plots() - each
took ownership of the whole figure and rebuilt it from scratch, so the picture
was a property of the last button pressed rather than of the data in memory.
That is why there was no way to see components again after switching views, and
no way to see copper and components at once.

The load-bearing part of that design was fig.clear(). plot_heatmap() created a
colorbar unconditionally, so the only way to stop colorbars accumulating was to
tear the figure down on every redraw. Tearing it down also destroyed the
AxesImage behind MainWindow.im, invalidated MainWindow.ax, and detached every
probe artist - which is where the H-5 and M-4 defects came from.

Here the axes are created once, artists are created once and then updated in
place, and the colorbar is retargeted rather than recreated. Nothing calls
fig.clear().

Placement note: this lives in its own module rather than in visualization.py
because visualization.py is deliberately Qt-free - tests/test_visualization.py
imports it without a QApplication, and pulling the Qt backend in there would
couple that test to a GUI stack it does not need.
"""

import logging

import matplotlib
import numpy as np
from PyQt6.QtCore import QTimer
from PyQt6.QtWidgets import QSizePolicy, QVBoxLayout, QWidget
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

import config

logger = logging.getLogger(__name__)


class BoardView(QWidget):
    """
    A single set of axes showing the board, composed of independent layers.

    Callers push data in (:meth:`set_temperature`, :meth:`set_copper`, ...) and
    mark what changed with :meth:`invalidate`. They never draw. Redraws are
    coalesced onto the next turn of the event loop, so a caller may invalidate
    the same layer any number of times and still pay for one repaint.
    """

    # --- Stacking order -----------------------------------------------------
    # Explicit constants rather than creation order: a translucent copper
    # overlay added after the component outlines would otherwise hide them.
    Z_BASE = 1
    Z_COPPER = 2
    Z_CONVECTION = 3
    Z_PARTS = 4
    Z_LABELS = 5
    Z_PROBES = 10

    #: Layers understood by :meth:`invalidate` and :meth:`set_layer_visible`.
    #: 'probes' holds artists created by the window and only registered here,
    #: so it takes part in visibility but has no renderer of its own.
    LAYERS = ("base", "copper", "convection", "components", "probes")

    #: Base fields are mutually exclusive: one full-size image, one colorbar.
    BASE_FIELDS = {
        "temperature": {"cmap": "inferno", "label": "Temperature [°C]"},
        "conductivity": {"cmap": "copper", "label": "Thermal Conductivity [W/mK]"},
        "power": {"cmap": "magma", "label": "Power density [W/m³]"},
        "none": {"cmap": None, "label": None},
    }

    #: Cells below this copper fraction are treated as bare laminate and left
    #: fully transparent, so the mask tints only what is actually copper.
    COPPER_MASK_THRESHOLD = 0.02
    COPPER_MASK_ALPHA = 0.45

    def __init__(self, parent=None, figsize=(8, 6)):
        super().__init__(parent)

        self.figure = Figure(figsize=figsize)
        self.canvas = FigureCanvas(self.figure)
        # Created once, for the lifetime of the widget. Never replaced.
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xlabel("X [mm]")
        self.ax.set_ylabel("Y [mm]")

        # The plot is the thing worth enlarging, so claim the spare space
        # rather than sharing it with the control strips below.
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)

        # --- Model: what there is to draw ---------------------------------
        self._temperature = None
        self._copper = None          # K matrix [W/mK]
        self._power = None           # Q matrix [W/m^3]
        self._convection = None      # h matrix [W/m^2K]
        self._components = []

        # --- View state: what to show -------------------------------------
        self._base_kind = "temperature"
        self._visible = {name: True for name in self.LAYERS}
        self._visible["copper"] = False
        self._visible["convection"] = False
        self._title = ""
        self._clim_locked = False
        self._clim_floor = None
        self._clim_running_max = None

        # --- Render state --------------------------------------------------
        self._artists = {}
        self._cbar = None
        self._dirty = set()
        self._redraw_pending = False

        self._copper_cmap = matplotlib.colormaps["copper"].copy()
        # Masked cells (bare laminate) must not tint the layer underneath.
        self._copper_cmap.set_bad(alpha=0.0)

    # ------------------------------------------------------------------
    # Data in. None of these draw anything.
    # ------------------------------------------------------------------

    def set_temperature(self, field):
        self._temperature = field
        if self._base_kind == "temperature":
            self.invalidate("base")

    def set_copper(self, k_matrix):
        """Conductivity field; feeds both the base view and the mask overlay."""
        self._copper = k_matrix
        self.invalidate("copper")
        if self._base_kind == "conductivity":
            self.invalidate("base")

    def set_power(self, q_matrix):
        self._power = q_matrix
        if self._base_kind == "power":
            self.invalidate("base")

    def set_convection(self, h_matrix):
        self._convection = h_matrix
        self.invalidate("convection")

    def set_components(self, components):
        self._components = list(components or ())
        self.invalidate("components")

    def set_title(self, text):
        self._title = text or ""
        self._schedule_redraw()

    def set_probe_artists(self, artists):
        """
        Registers externally created probe artists for visibility control.

        The window creates them - it knows the positions and the readings - but
        the view decides whether they are shown, so every layer is toggled the
        same way.
        """
        self._artists["probes"] = list(artists or ())
        self._schedule_redraw()

    # ------------------------------------------------------------------
    # View state
    # ------------------------------------------------------------------

    def set_base_kind(self, kind):
        if kind not in self.BASE_FIELDS:
            raise ValueError(
                f"Unknown base field '{kind}'. "
                f"Known: {', '.join(sorted(self.BASE_FIELDS))}."
            )
        if kind == self._base_kind:
            return
        self._base_kind = kind
        self._clim_running_max = None
        self.invalidate("base")

    def set_layer_visible(self, name, visible):
        if name not in self._visible:
            raise ValueError(f"Unknown layer '{name}'. Known: {', '.join(self.LAYERS)}.")
        if self._visible[name] == bool(visible):
            return
        self._visible[name] = bool(visible)
        # Visibility is cheap: no data rebuild, just a repaint.
        self._schedule_redraw()

    def is_layer_visible(self, name):
        return self._visible[name]

    def set_clim_lock(self, locked, floor=None):
        """
        Pins the colour scale so a transient can be read.

        Rescaling to each frame's own min/max makes the animation meaningless:
        the same colour is a different temperature on every frame. Locked, the
        bottom is held at `floor` (ambient) and the top only ever grows.
        """
        self._clim_locked = bool(locked)
        self._clim_floor = None if floor is None else float(floor)
        if not locked:
            self._clim_running_max = None
        self.invalidate("base")

    def reset_clim(self):
        """Forgets the running maximum, e.g. when a new run starts."""
        self._clim_running_max = None
        self.invalidate("base")

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def base_array(self, kind=None):
        """The array a given base kind would display, or None."""
        return {
            "temperature": self._temperature,
            "conductivity": self._copper,
            "power": self._power,
            "none": None,
        }[kind or self._base_kind]

    def field(self):
        """
        The array currently displayed as the base layer.

        Single source of truth for anything that reads values off the screen -
        the probes in particular. Keeping a second copy in the window is what
        lets a probe report a number that is not the one being shown.
        """
        return self.base_array()

    def field_kind(self):
        return self._base_kind

    def has_field(self):
        return self.base_array() is not None

    def has_copper(self):
        return self._copper is not None and float(np.ptp(self._copper)) > 0.0

    def has_convection_zones(self):
        return self._convection is not None and float(np.ptp(self._convection)) > 0.0

    # ------------------------------------------------------------------
    # Redraw lifecycle
    # ------------------------------------------------------------------

    def invalidate(self, *layers):
        """
        Marks layers as needing their data rebuilt, then schedules one repaint.

        Called with no arguments, every layer is rebuilt.
        """
        self._dirty.update(layers or self.LAYERS)
        self._schedule_redraw()

    def request_draw(self):
        """Repaint without rebuilding anything (used after external artists)."""
        self._schedule_redraw()

    def flush_pending(self):
        """
        Applies any deferred rebuild right now.

        Needed before anything that reads the figure synchronously - savefig()
        renders from the artists, so a pending rebuild would be exported stale.
        """
        if self._redraw_pending or self._dirty:
            self._flush()

    def _schedule_redraw(self):
        # Coalesces a burst of changes into a single repaint. This is what
        # makes the blockSignals guard around add_component_row unnecessary:
        # six cellChanged signals now cost one redraw, not seven.
        if self._redraw_pending:
            return
        self._redraw_pending = True
        QTimer.singleShot(0, self._flush)

    def _flush(self):
        self._redraw_pending = False
        try:
            dirty, self._dirty = self._dirty, set()
            for name in dirty:
                self._rebuild(name)
            self._apply_visibility()
            self.ax.set_title(self._title)
            self.canvas.draw_idle()
        except RuntimeError:
            # Widget torn down between scheduling and firing.
            pass

    def _rebuild(self, name):
        if name == "base":
            self._rebuild_base()
        elif name == "copper":
            self._rebuild_copper()
        elif name == "convection":
            self._rebuild_convection()
        elif name == "components":
            self._rebuild_components()
        # 'probes' has no renderer: the window supplies those artists.

    def _base_is_showable(self):
        return (self._visible.get("base", True)
                and self._base_kind != "none"
                and self.base_array() is not None)

    def _apply_visibility(self):
        effective = dict(self._visible)
        effective["base"] = self._base_is_showable()
        # An overlay with nothing to draw stays hidden regardless of intent.
        if not self.has_copper():
            effective["copper"] = False
        if not self.has_convection_zones():
            effective["convection"] = False

        for name, visible in effective.items():
            for artist in self._artists_of(name):
                artist.set_visible(visible)
        if self._cbar is not None:
            self._cbar.ax.set_visible(effective["base"])

    def _artists_of(self, name):
        entry = self._artists.get(name)
        if entry is None:
            return ()
        return entry if isinstance(entry, (list, tuple)) else (entry,)

    # ------------------------------------------------------------------
    # Layer renderers
    # ------------------------------------------------------------------

    def _extent(self):
        """
        Shared by every image layer, so overlays cannot misregister.

        Axes are in CAD millimetres, not board-relative ones: the component
        coordinates in the CSV and the pads in the Gerber are both quoted in
        the CAD frame, and a reader comparing the plot to the schematic needs
        the same numbers.
        """
        x0, y0, x1, y1 = config.board_extent()
        return [x0, x1, y0, y1]

    def _rebuild_base(self):
        data = self.base_array()
        if data is None or self._base_kind == "none":
            return
        spec = self.BASE_FIELDS[self._base_kind]

        im = self._artists.get("base")
        if im is None:
            im = self.ax.imshow(
                data, origin="lower", extent=self._extent(), zorder=self.Z_BASE,
            )
            self._artists["base"] = im
        else:
            im.set_data(data)
            im.set_extent(self._extent())
        im.set_cmap(spec["cmap"])
        im.set_clim(*self._clim_for(data))

        # Created once; retargeted afterwards. Creating a second colorbar is
        # what forced fig.clear() in the old design.
        if self._cbar is None:
            self._cbar = self.figure.colorbar(im, ax=self.ax)
        else:
            self._cbar.update_normal(im)
        self._cbar.set_label(spec["label"])

    def _clim_for(self, data):
        vmin = float(np.nanmin(data))
        vmax = float(np.nanmax(data))
        if not np.isfinite(vmin) or not np.isfinite(vmax):
            vmin, vmax = 0.0, 1.0

        if self._clim_locked and self._base_kind == "temperature":
            if self._clim_floor is not None:
                vmin = self._clim_floor
            self._clim_running_max = (
                vmax if self._clim_running_max is None
                else max(self._clim_running_max, vmax)
            )
            vmax = self._clim_running_max

        if vmax <= vmin:
            # A uniform field - the ambient reset before a run - has no range.
            # The pad has to be visible, not epsilon: nudging by 1e-9 makes
            # Matplotlib fall back to offset notation and label the bar
            # "1e-9+2.5e1" instead of a readable span around 25 °C.
            pad = max(abs(vmin) * 0.02, 0.5)
            vmin, vmax = vmin - pad, vmin + pad
        return vmin, vmax

    def _copper_coverage(self):
        """
        Copper fraction per cell, derived from the conductivity field.

        Inverts the model from FIX C-1: k_cell = k_base + coverage * k_delta,
        so coverage falls straight out. Because the Gerber loader area-averages
        (FIX M-5), this is a genuine fractional value - dense pours come out
        more opaque than thin traces, and the transparency carries information
        rather than just decorating.
        """
        if self._copper is None:
            return None
        k_base = config.calculate_k_eff(
            layers=config.BOARD_LAYERS, copper_oz=config.COPPER_OZ,
            substrate_k=config.K_FR4,
        )
        k_delta = (config.K_CU * config.COPPER_OZ * 35e-6) / config.d
        if k_delta <= 0:
            return None
        coverage = np.clip((self._copper - k_base) / k_delta, 0.0, 1.0)
        return np.ma.masked_where(coverage < self.COPPER_MASK_THRESHOLD, coverage)

    def _rebuild_copper(self):
        coverage = self._copper_coverage()
        if coverage is None:
            return
        im = self._artists.get("copper")
        if im is None:
            im = self.ax.imshow(
                coverage, cmap=self._copper_cmap, alpha=self.COPPER_MASK_ALPHA,
                vmin=0.0, vmax=1.0, origin="lower", extent=self._extent(),
                zorder=self.Z_COPPER,
            )
            self._artists["copper"] = im
        else:
            im.set_data(coverage)
            im.set_extent(self._extent())

    def _rebuild_convection(self):
        for artist in self._artists_of("convection"):
            try:
                artist.remove()
            except (ValueError, NotImplementedError, AttributeError):
                pass
        self._artists["convection"] = []

        if not self.has_convection_zones():
            return
        h = self._convection
        # Outlines, not a fill: heatsink zones sit on top of the components
        # they cool, and a filled patch would bury them.
        level = float(h.min()) + (float(h.max()) - float(h.min())) * 0.5
        contours = self.ax.contour(
            h, levels=[level], colors="cyan", linestyles=":", linewidths=1.3,
            origin="lower", extent=self._extent(), zorder=self.Z_CONVECTION,
        )
        self._artists["convection"] = contours

    def _rebuild_components(self):
        for artist in self._artists_of("components"):
            try:
                artist.remove()
            except (ValueError, NotImplementedError):
                pass
        artists = []

        board_x0, board_y0, board_x1, board_y1 = config.board_extent()
        for comp in self._components:
            width = comp.get("Width_mm") or 0.0
            length = comp.get("Length_mm") or 0.0
            if width <= 0 or length <= 0:
                continue
            x0 = comp["Center_X_mm"] - width / 2.0
            y0 = comp["Center_Y_mm"] - length / 2.0
            on_board = (x0 + width > board_x0 and x0 < board_x1
                        and y0 + length > board_y0 and y0 < board_y1)
            colour = "deepskyblue" if on_board else "red"
            rect = Rectangle(
                (x0, y0), width, length, fill=False, edgecolor=colour,
                linewidth=1.2, linestyle="-" if on_board else "--",
                zorder=self.Z_PARTS,
            )
            self.ax.add_patch(rect)
            artists.append(rect)
            artists.append(self.ax.annotate(
                f"{comp['Designator']}\n{comp.get('Power_Watts', 0.0):.2f} W",
                xy=(comp["Center_X_mm"], comp["Center_Y_mm"]),
                ha="center", va="center", fontsize=7, color=colour,
                zorder=self.Z_LABELS, clip_on=True,
            ))

        self._artists["components"] = artists
