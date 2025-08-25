# gaussbeam2d_fit.py
# PyQt5 port of the MATLAB "Gaussbeam2DFit" GUI and fitting workflow.

import sys
import math
import numpy as np
from dataclasses import dataclass, asdict
from typing import List, Tuple, Optional
from matplotlib.widgets import RectangleSelector

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QCheckBox, QComboBox, QFileDialog, QMessageBox,
    QInputDialog, QTableWidget, QTableWidgetItem, QDialog, QFormLayout,
    QGroupBox, QSplitter, QLineEdit, QFrame, QScrollArea
)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
import matplotlib.pyplot as plt

from PIL import Image
from scipy.optimize import fmin, minimize


# ---------------------------- Math helpers ---------------------------- #

def rotated_gaussian(params, x, y, offset_on):
    """
    params = [w, wbx, wby, Imax, cx, cy] (+ [offst] if offset_on)
    w is rotation in radians
    wbx, wby are 4*sigma (D4σ) widths in pixels
    """
    w, wbx, wby, Imax, cx, cy = params[:6]
    offst = params[6] if (offset_on and len(params) >= 7) else 0.0

    cosw = np.cos(w)
    sinw = np.sin(w)
    xp = (cosw * (x - cx) + sinw * (y - cy))
    yp = (-sinw * (x - cx) + cosw * (y - cy))

    # MATLAB code uses exp(-2 * (xp/wbx)^2 - 2 * (yp/wby)^2)
    f = Imax * np.exp(-2.0 * (xp ** 2) / (wbx ** 2) - 2.0 * (yp ** 2) / (wby ** 2))
    return f - offst


def sse_for_image(params, x, y, img, offset_on):
    f = rotated_gaussian(params, x, y, offset_on)
    diff = f - img
    return np.sum(diff * diff)


def normalized_chi(f, img):
    num = np.sum((f - img) ** 2)
    den = np.sum(f ** 2) + 1e-12
    return float(num / den)


def wrap_degrees(d):
    # Map angle to (-180, 180]
    if np.isfinite(d):
        while d > 180:
            d -= 360
        while d <= -180:
            d += 360
    return d


def fwhm_1d(profile):
    """Return index width at half max around the maximum."""
    profile = np.asarray(profile, dtype=float)
    idx_max = int(np.argmax(profile))
    half = 0.5 * profile[idx_max]

    # left crossing
    i_left = idx_max
    while i_left > 0 and profile[i_left] > half:
        i_left -= 1
    # right crossing
    i_right = idx_max
    n = len(profile)
    while i_right < n - 1 and profile[i_right] > half:
        i_right += 1
    return max(1, i_right - i_left)  # number of samples between crossings


def bilinear_sample(img, x, y):
    """Bilinear sampling of 2D array img at fractional coords (x,y)."""
    h, w = img.shape
    x0 = np.floor(x).astype(int)
    y0 = np.floor(y).astype(int)
    x1 = x0 + 1
    y1 = y0 + 1

    x0 = np.clip(x0, 0, w - 1)
    x1 = np.clip(x1, 0, w - 1)
    y0 = np.clip(y0, 0, h - 1)
    y1 = np.clip(y1, 0, h - 1)

    Ia = img[y0, x0]
    Ib = img[y1, x0]
    Ic = img[y0, x1]
    Id = img[y1, x1]

    wa = (x1 - x) * (y1 - y)
    wb = (x1 - x) * (y - y0)
    wc = (x - x0) * (y1 - y)
    wd = (x - x0) * (y - y0)

    return wa * Ia + wb * Ib + wc * Ic + wd * Id


# ---------------------------- Data structures ---------------------------- #

@dataclass
class FitResult:
    chi: float
    x: np.ndarray
    y: np.ndarray
    f_fit: np.ndarray
    wbx: float
    wby: float
    Imax: float
    pic: np.ndarray
    offst: float
    angle_deg: float
    cx: float
    cy: float
    angle_deg_raw: float  # unwrapped before [-180,180] mapping


@dataclass
class ImageEntry:
    fit: FitResult
    picture: np.ndarray
    dist_um: float  # distance to starting point (microns)


# ---------------------------- UI: Matplotlib canvases ---------------------------- #

class MplCanvas(FigureCanvas):
    def __init__(self, width=4, height=3, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)
        self.ax1 = self.fig.add_subplot(1, 2, 1)
        self.ax2 = self.fig.add_subplot(1, 2, 2)
        self.fig.tight_layout()


class MplSingleCanvas(FigureCanvas):
    def __init__(self, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)
        self.ax = self.fig.add_subplot(1, 1, 1)
        self.fig.tight_layout()


# ---------------------------- Crop dialog ---------------------------- #

class CropDialog(QDialog):
    def __init__(self, img_np, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select a square (rectangle will be snapped to square)")
        self.img_np = img_np
        self.crop_rect = None

        self.canvas = MplSingleCanvas(width=6, height=6, dpi=100)
        self.canvas.ax.imshow(img_np, cmap='gray')
        self.canvas.ax.set_title("Drag to select a region")

        # FIX: remove drawtype argument
        self.selector = RectangleSelector(
            self.canvas.ax,
            self.on_select,
            useblit=True,
            interactive=True
        )

        btn_box = QHBoxLayout()
        ok_btn = QPushButton("OK")
        cancel_btn = QPushButton("Cancel")
        ok_btn.clicked.connect(self.accept)
        cancel_btn.clicked.connect(self.reject)
        btn_box.addStretch(1)
        btn_box.addWidget(ok_btn)
        btn_box.addWidget(cancel_btn)

        layout = QVBoxLayout(self)
        layout.addWidget(self.canvas)
        # Toolbar for save/zoom in crop dialog
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addLayout(btn_box)
        self.setLayout(layout)
        self.resize(700, 700)

    def on_select(self, eclick, erelease):
        x0, y0 = int(eclick.xdata), int(eclick.ydata)
        x1, y1 = int(erelease.xdata), int(erelease.ydata)
        x_min, x_max = sorted([x0, x1])
        y_min, y_max = sorted([y0, y1])
        w = x_max - x_min
        h = y_max - y_min
        if w <= 0 or h <= 0:
            return
        side = min(w, h)
        self.crop_rect = (x_min, y_min, side, side)

    def get_cropped(self) -> Optional[np.ndarray]:
        if not self.crop_rect:
            return None
        x, y, s, _ = self.crop_rect
        y2 = y + s
        x2 = x + s
        return self.img_np[y:y2, x:x2].copy()


# ---------------------------- Main Window ---------------------------- #

class GaussBeamApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Gaussbeam 2D Fit (PyQt)")
        self.resize(1200, 720)

        # State
        self.pixel_size_um = 3.75  # default pixelsize in microns
        self.offset_on = False
        self.wR_rad = 0.0  # angle guess in radians
        self.entries: List[ImageEntry] = []  # list of images + fits
        # self.cbar1 = None
        # self.cbar2 = None

        # Central UI
        central = QWidget()
        self.setCentralWidget(central)
        outer = QHBoxLayout(central)

        # Left: plots
        self.canvas = MplCanvas(width=7, height=5, dpi=100)

        # Right: controls + table
        right = QVBoxLayout()
        btn_row1 = QHBoxLayout()
        self.btn_intensity = QPushButton("Intensity Profile")
        self.btn_add = QPushButton("Add Image")
        self.btn_save = QPushButton("Save Results")
        btn_row1.addWidget(self.btn_intensity)
        btn_row1.addWidget(self.btn_add)
        btn_row1.addWidget(self.btn_save)

        # Offset checkbox
        self.chk_offset = QCheckBox("Offset")
        self.chk_offset.setChecked(False)

        # Pixel size selector
        pix_row = QHBoxLayout()
        pix_row.addWidget(QLabel("Pixelsize [µm]"))
        self.cbo_pix = QComboBox()
        self.cbo_pix.addItems(["3.75", "3.10", "3.69", "4.54", "5.55", "5.86", "6.45"])
        self.cbo_pix.setCurrentText("3.75")
        pix_row.addWidget(self.cbo_pix)
        pix_row.addStretch(1)

        # Image selector row  
        sel_row = QHBoxLayout()
        sel_row.addWidget(QLabel("Show image #"))
        self.cbo_image = QComboBox()
        self.cbo_image.setMinimumWidth(90)
        self.cbo_image.currentIndexChanged.connect(self.on_image_selected)
        self.btn_prev = QPushButton("<")
        self.btn_next = QPushButton(">")
        self.btn_prev.clicked.connect(self.on_prev_image)
        self.btn_next.clicked.connect(self.on_next_image)

        self.btn_delete = QPushButton("Delete")          
        self.btn_delete.clicked.connect(self.on_delete_image)  

        sel_row.addWidget(self.cbo_image)
        sel_row.addWidget(self.btn_prev)
        sel_row.addWidget(self.btn_next)
        sel_row.addWidget(self.btn_delete) 
        sel_row.addStretch(1)
        # Parameter table
        self.tbl = QTableWidget(6, 1)
        self.tbl.setHorizontalHeaderLabels(["values"])
        self.tbl.setVerticalHeaderLabels([
            "waist x [µm]", "waist y [µm]", "angle(x) [°]",
            "maximum intensity", "normalized chi-square", "offset"
        ])
        self.tbl.horizontalHeader().setStretchLastSection(True)

        # Below: "Results" actions like MATLAB's second window
        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setFrameShadow(QFrame.Sunken)

        res_row = QHBoxLayout()
        self.btn_add_next = QPushButton("Add Next Image")
        self.btn_profile = QPushButton("Create Profile")
        self.btn_new = QPushButton("New Measurement")
        self.btn_save_all = QPushButton("Save All")
        res_row.addWidget(self.btn_add_next)
        res_row.addWidget(self.btn_profile)
        res_row.addWidget(self.btn_new)
        res_row.addWidget(self.btn_save_all)
        res_row.addStretch(1)

        # Wire up
        self.btn_add.clicked.connect(self.on_add_image_initial)
        self.btn_save.clicked.connect(self.on_save_single)
        self.chk_offset.stateChanged.connect(self.on_toggle_offset)
        self.cbo_pix.currentTextChanged.connect(self.on_pixel_size_changed)

        self.btn_add_next.clicked.connect(self.on_add_image_next)
        self.btn_profile.clicked.connect(self.on_create_profile)
        self.btn_new.clicked.connect(self.on_new_measurement)
        self.btn_save_all.clicked.connect(self.on_save_all)
        self.btn_intensity.clicked.connect(self.on_intensity_profile)

        # Layout pack
        left_box = QVBoxLayout()
        left_box.addWidget(self.canvas)
        # Navigation toolbar for Save/Zoom/Pan
        self.toolbar_main = NavigationToolbar(self.canvas, self)
        left_box.addWidget(self.toolbar_main)

        right.addLayout(btn_row1)
        right.addWidget(self.chk_offset)
        right.addLayout(pix_row)
        right.addLayout(sel_row) 
        right.addWidget(self.tbl)
        right.addWidget(sep)
        right.addLayout(res_row)
        right.addStretch(1)

        outer.addLayout(left_box, 3)
        outer.addLayout(right, 2)

        # Prompt to load first image
        QtCore.QTimer.singleShot(200, self.on_add_image_initial)

    # Function added to refresh the picture

    def on_delete_image(self):
        if not self.entries:
            return
        idx = self.current_index
        # Confirm
        reply = QMessageBox.question(
            self,
            "Delete image",
            f"Delete image #{idx + 1}?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )
        if reply != QMessageBox.Yes:
            return

        # Remove it
        del self.entries[idx]

        if not self.entries:
            # Nothing left: clear UI
            self.current_index = 0
            self.canvas.ax1.clear()
            self.canvas.ax2.clear()
            self.canvas.draw_idle()
            self.tbl.clearContents()
            self.refresh_image_selector()
            return

        # Clamp index to last element if needed
        if self.current_index >= len(self.entries):
            self.current_index = len(self.entries) - 1

        # Show current image and refresh selector
        self.show_image_by_index(self.current_index)
        

    def refresh_image_selector(self):
        self.cbo_image.blockSignals(True)
        self.cbo_image.clear()
        if self.entries:
            self.cbo_image.addItems([f"{i+1}" for i in range(len(self.entries))])
            self.cbo_image.setCurrentIndex(self.current_index)
        self.cbo_image.blockSignals(False)

        # enable/disable prev/next
        enabled = len(self.entries) > 0
        self.btn_prev.setEnabled(enabled and self.current_index > 0)
        self.btn_next.setEnabled(enabled and self.current_index < len(self.entries) - 1)
        self.btn_delete.setEnabled(enabled)   

    def show_image_by_index(self, idx: int):
        if not self.entries:
            return
        idx = max(0, min(idx, len(self.entries) - 1))
        self.current_index = idx
        fr = self.entries[idx].fit
        self.update_plots_and_table(fr)
        self.refresh_image_selector()

    def on_image_selected(self, idx: int):
        if not self.entries:
            return
        self.show_image_by_index(idx)

    def on_prev_image(self):
        if self.entries and self.current_index > 0:
            self.show_image_by_index(self.current_index - 1)

    def on_next_image(self):
        if self.entries and self.current_index < len(self.entries) - 1:
            self.show_image_by_index(self.current_index + 1)
    # -------------------- Image I/O and cropping -------------------- #

    def read_image_and_crop(self) -> Optional[np.ndarray]:
        path, _ = QFileDialog.getOpenFileName(self, "Open image", "", "Images (*.png *.jpg *.jpeg *.tif *.tiff *.bmp)")
        if not path:
            return None
        # Read grayscale
        try:
            im = Image.open(path)
            if im.mode not in ("L", "I;16", "F"):
                im = im.convert("L")
            img = np.array(im).astype(np.float64)
            # Normalize to 0..1
            if img.max() > 0:
                img = img / img.max()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to read image:\n{e}")
            return None

        # Crop dialog
        QMessageBox.information(self, "Read in image", "Select a square. Rectangle will be cut off to square.")
        dlg = CropDialog(img, self)
        if dlg.exec_() != QDialog.Accepted:
            return None
        cropped = dlg.get_cropped()
        if cropped is None or cropped.size == 0:
            QMessageBox.warning(self, "No selection", "No crop was selected.")
            return None

        # Ensure square (safety)
        h, w = cropped.shape
        s = min(h, w)
        cropped = cropped[:s, :s]
        return cropped

    # -------------------- Fitting pipeline -------------------- #

    def create_fit_for_image(self, picture: np.ndarray) -> FitResult:
        # square to min dimension already handled; grid:
        n = min(picture.shape[0], picture.shape[1])
        x = np.arange(n)
        y = np.arange(n)
        X, Y = np.meshgrid(x, y)

        pic = picture[:n, :n]

        params0 = self.find_start_parameters(pic, self.wR_rad)
        if not self.offset_on:
            p0 = [params0['w'], params0['wbx'], params0['wby'], params0['Imax'], params0['cx'], params0['cy']]
        else:
            p0 = [params0['w'], params0['wbx'], params0['wby'], params0['Imax'], params0['cx'], params0['cy'], params0['offst']]

        # Use Nelder-Mead like MATLAB fminsearch
        def obj(p):
            return sse_for_image(p, X, Y, pic, self.offset_on)

        p_opt = fmin(obj, p0, maxiter=100000, maxfun=100000, disp=False)

        # Build result and chi
        f = rotated_gaussian(p_opt, X, Y, self.offset_on)
        chi = normalized_chi(f, pic)

        w_rad = float(p_opt[0])
        w_deg_raw = w_rad * 180.0 / math.pi
        w_deg = wrap_degrees(w_deg_raw)
        self.wR_rad = w_rad  # update guess for next image

        wbx, wby, Imax, cx, cy = float(p_opt[1]), float(p_opt[2]), float(p_opt[3]), float(p_opt[4]), float(p_opt[5])
        offst = float(p_opt[6]) if self.offset_on and len(p_opt) >= 7 else 0.0

        return FitResult(
            chi=chi, x=X, y=Y, f_fit=f, wbx=wbx, wby=wby, Imax=Imax,
            pic=pic, offst=offst, angle_deg=w_deg, cx=cx, cy=cy, angle_deg_raw=w_deg_raw
        )

    def find_start_parameters(self, img: np.ndarray, wR_rad: float):
        # sum over two dimensions
        prof_x = img.sum(axis=0)  # along rows -> columns (x)
        prof_y = img.sum(axis=1)  # along columns -> rows (y)
        cx = int(np.argmax(prof_x))
        cy = int(np.argmax(prof_y))
        Imax = float(img[cy, cx])
        offst = float(img[0, 0])

        # FWHM around center (scan along row/col that pass through cy/cx)
        widthx = fwhm_1d(img[cy, :])
        widthy = fwhm_1d(img[:, cx])

        # MATLAB initializes w from previous result
        w = float(wR_rad)

        return dict(cx=float(cx), cy=float(cy), Imax=Imax,
                    wbx=float(widthx), wby=float(widthy),
                    w=w, offst=offst)

    # -------------------- UI Updates -------------------- #

    def update_plots_and_table(self, fr: FitResult):
        ps = self.pixel_size_um

        # --- Clear previous content ---
        self.canvas.ax1.clear()
        self.canvas.ax2.clear()

        extent = [0, fr.x.shape[1] * ps, 0, fr.y.shape[0] * ps]

        # --- Plot image ---
        im1 = self.canvas.ax1.imshow(fr.pic, extent=extent, origin='lower', cmap='gray')
        self.canvas.ax1.set_title("image")
        self.canvas.ax1.set_xlabel("x [µm]")
        self.canvas.ax1.set_ylabel("y [µm]")

        # --- Plot fit ---
        im2 = self.canvas.ax2.imshow(fr.f_fit, extent=extent, origin='lower', cmap='inferno')
        self.canvas.ax2.set_title("fit")
        self.canvas.ax2.set_xlabel("x [µm]")
        self.canvas.ax2.set_ylabel("y [µm]")

        # Redraw canvas
        self.canvas.draw_idle()

        extent = [0, fr.x.shape[1] * ps, 0, fr.y.shape[0] * ps]


        # Table
        rows = [
            fr.wbx * ps,
            fr.wby * ps,
            fr.angle_deg,
            fr.Imax,
            fr.chi,
            fr.offst
        ]
        for i, val in enumerate(rows):
            item = QTableWidgetItem(f"{val:.6g}")
            self.tbl.setItem(i, 0, item)

    # -------------------- Buttons -------------------- #

    def on_add_image_initial(self):
        pic = self.read_image_and_crop()
        if pic is None:
            return
        fr = self.create_fit_for_image(pic)
        entry = ImageEntry(fit=fr, picture=pic, dist_um=0.0)
        self.entries = [entry]
        self.current_index = 0
        self.update_plots_and_table(fr)
        self.refresh_image_selector()                 # <-- ensure selector shows "1"
        QMessageBox.information(self, "Loaded", "First image loaded and fitted.")

    def on_add_image_next(self):
        if not self.entries:
            self.on_add_image_initial()
            return
        pic = self.read_image_and_crop()
        if pic is None:
            return
        fr = self.create_fit_for_image(pic)

        dist_mm, ok = QInputDialog.getDouble(self, "Distance", "Enter distance to starting point [mm]:", 0.0, -1e9, 1e9, 3)
        if not ok:
            return
        dist_um = float(dist_mm) * 1000.0
        self.entries.append(ImageEntry(fit=fr, picture=pic, dist_um=dist_um))
        self.current_index = len(self.entries) - 1     # jump to the new image
        self.update_plots_and_table(fr)
        self.refresh_image_selector()
        QMessageBox.information(self, "Added", f"Image #{len(self.entries)} added (distance {dist_mm} mm).")
    def on_save_single(self):
        if not self.entries:
            return
        fr = self.entries[0].fit
        ps = self.pixel_size_um
        # MATLAB "beamFitFile1.txt" (simple)
        path, _ = QFileDialog.getSaveFileName(self, "Save Results", "beamFitFile1.txt", "Text (*.txt)")
        if not path:
            return
        header = "imagenumber normalized_chi-square waist_x waist_y I_max Offset Angle\n"
        vals = [1, fr.chi, fr.wbx * ps, fr.wby * ps, fr.Imax, fr.offst, fr.angle_deg]
        with open(path, "w", encoding="utf-8") as f:
            f.write(header)
            f.write("{:6.0f} {:6.4f} {:6.4f} {:6.4f} {:6.4f} {:6.4f} {:6.4f}\n".format(*vals))
        QMessageBox.information(self, "Saved", f"Saved to {path}")

    def on_save_all(self):
        if not self.entries:
            return
        ps = self.pixel_size_um
        path, _ = QFileDialog.getSaveFileName(self, "Save All", "beamFitFileAll.txt", "Text (*.txt)")
        if not path:
            return
        header = "imagenumber normalized_chi-square waist_x waist_y I_max Offset Angle\n"
        with open(path, "w", encoding="utf-8") as f:
            f.write(header)
            for i, e in enumerate(self.entries, start=1):
                fr = e.fit
                vals = [i, fr.chi, fr.wbx * ps, fr.wby * ps, fr.Imax, fr.offst, fr.angle_deg]
                f.write("{:6.0f} {:6.4f} {:6.4f} {:6.4f} {:6.4f} {:6.4f} {:6.4f}\n".format(*vals))
        QMessageBox.information(self, "Saved", f"Saved all to {path}")

    def on_toggle_offset(self, state):
        self.offset_on = (state == QtCore.Qt.Checked)
        if self.entries:
            # Refit current picture with new offset setting
            pic = self.entries[-1].picture
            fr = self.create_fit_for_image(pic)
            self.entries[-1] = ImageEntry(fit=fr, picture=pic, dist_um=self.entries[-1].dist_um)
            self.update_plots_and_table(fr)

    def on_pixel_size_changed(self, text):
        try:
            self.pixel_size_um = float(text)
        except:
            return
        if self.entries:
            self.update_plots_and_table(self.entries[-1].fit)

    def on_new_measurement(self):
        self.entries = []
        self.wR_rad = 0.0
        self.current_index = 0 
        self.canvas.ax1.clear(); self.canvas.ax2.clear(); self.canvas.draw_idle()
        self.tbl.clearContents()
        self.refresh_image_selector()
        self.on_add_image_initial()

    # -------------------- Intensity profiles -------------------- #

    def on_intensity_profile(self):
        if not self.entries:
            return
        # always use the image currently selected in the “Show image #” selector
        self.show_intensity_profile_for_index(self.current_index)


    def show_intensity_profile_for_index(self, idx: int):
        e = self.entries[idx]
        fr = e.fit

        # Build two orthogonal lines through (cx,cy) at angle w and w+90°
        # Sample along the line until reaching image boundary.
        h, w = fr.pic.shape
        cx, cy = fr.cx, fr.cy
        ang_deg = fr.angle_deg
        ang_rad = ang_deg * math.pi / 180.0

        def sample_line(theta):
            # direction vector
            dx = math.cos(theta)
            dy = math.sin(theta)
            # step out to bounds
            # estimate max t
            ts = []
            for sign in (-1, +1):
                t = 0.0
                x = cx
                y = cy
                # large enough steps; stop when beyond bounds
                while 1 <= x < w - 2 and 1 <= y < h - 2:
                    x += sign * dx
                    y += sign * dy
                    t += 1.0
                ts.append(t)
            tmax = int(max(ts))
            ts_full = np.linspace(-tmax, tmax, 2 * tmax + 1)
            xs = cx + ts_full * dx
            ys = cy + ts_full * dy
            data = bilinear_sample(fr.pic, xs, ys)
            model = bilinear_sample(fr.f_fit, xs, ys)
            return ts_full, data, model

        t1, d1, m1 = sample_line(ang_rad)  # x-waist direction
        t2, d2, m2 = sample_line(ang_rad + math.pi / 2.0)  # y-waist

        dlg = QDialog(self)
        dlg.setWindowTitle(f"Intensity measurement and fit (image #{idx + 1})")
        v = QVBoxLayout(dlg)

        canvas = FigureCanvas(Figure(dpi=100, layout='constrained'))
        ax1 = canvas.figure.add_subplot(2, 1, 1)
        ax2 = canvas.figure.add_subplot(2, 1, 2)

        ax1.plot(t1, d1, "o", label="data")
        ax1.plot(t1, m1, label="fit")
        ax1.set_title("x-waist direction")
        # ax1.set_xlabel("samples")
        ax1.set_ylabel("Intensity")
        ax1.legend()

        ax2.plot(t2, d2, "o", label="data")
        ax2.plot(t2, m2, label="fit")
        ax2.set_title("y-waist direction")
        ax2.set_xlabel("samples")
        ax2.set_ylabel("Intensity")
        ax2.legend()

        v.addWidget(canvas)
        toolbar = NavigationToolbar(canvas, dlg)
        v.addWidget(toolbar)
        dlg.resize(800, 600)
        dlg.exec_()

    # -------------------- Propagation profile (w(z)) -------------------- #

    def on_create_profile(self):
        if not self.entries:
            return

        # Ask wavelength (nm)
        wl_str, ok = QInputDialog.getText(self, "Wavelength", "Enter the wavelength [nm]:", text="780")
        if not ok:
            return
        try:
            wavelength_um = float(wl_str) / 1000.0  # convert nm -> µm
        except:
            QMessageBox.warning(self, "Invalid input", "Could not parse wavelength.")
            return

        # Build arrays for x and y waists (in µm) and z positions (in µm)
        ps = self.pixel_size_um
        widths_x_pix = []
        widths_y_pix = []
        z_list_um = [0.0]

        for i, e in enumerate(self.entries):
            fr = e.fit
            widths_x_pix.append(fr.wbx)
            widths_y_pix.append(fr.wby)
            if i == 0:
                continue
            z_list_um.append(float(e.dist_um))

        widths_x = np.array(widths_x_pix) * ps
        widths_y = np.array(widths_y_pix) * ps
        z_um = np.array(z_list_um)

        # Fit function: w(z) = w0 * sqrt(1 + M^4 * ((z - z0) / zR)^2)
        # where zR = pi * w0^2 / lambda
        def w_model(z, w0, z0, M, wl_um):
            zR = math.pi * (w0 ** 2) / wl_um
            return w0 * np.sqrt(1.0 + (M ** 4) * ((z - z0) / zR) ** 2)

        # Objective for a given dimension (x or y)
        def fit_profile(meas_w_um, z_um):
            if len(meas_w_um) < 2:
                # not enough points: return dummy
                return dict(w0=meas_w_um.min(), z0=0.0, M=1.0, chi=0.0, zR=math.pi * (meas_w_um.min() ** 2) / wavelength_um)

            # Initial guesses
            w0_init = float(np.min(meas_w_um))
            z0_init = float(z_um[1] if len(z_um) > 1 else 0.0)
            M_init = 5.0

            def obj(p):
                w0, z0, M = p
                pred = w_model(z_um, w0, z0, M, wavelength_um)
                return np.sum((pred - meas_w_um) ** 2)

            bounds = [(1e-6, None), (None, None), (1.0, None)]  # w0>0, M>=1
            res = minimize(obj, x0=[w0_init, z0_init, M_init], method='SLSQP', bounds=bounds,
                           options=dict(maxiter=100000))

            w0, z0, M = res.x
            pred_points = w_model(z_um, w0, z0, M, wavelength_um)
            chi = float(np.sum((pred_points - meas_w_um) ** 2) / (np.sum(pred_points ** 2) + 1e-12))
            zR = math.pi * (w0 ** 2) / wavelength_um
            return dict(w0=float(w0), z0=float(z0), M=float(M), chi=chi, zR=float(zR))

        fit_x = fit_profile(widths_x, z_um)
        fit_y = fit_profile(widths_y, z_um)

        # Plot profile
        dlg = QDialog(self)
        dlg.setWindowTitle("Profiles (w(z))")
        layout = QVBoxLayout(dlg)

        canvas = FigureCanvas(Figure(figsize=(6, 7), dpi=100, layout='constrained'))
        ax1 = canvas.figure.add_subplot(2, 1, 1)
        ax2 = canvas.figure.add_subplot(2, 1, 2)

        z_line_um = np.linspace(0.0, z_um.max() if z_um.size > 0 else 1.0, 1000)
        ax1.plot(z_line_um / 1000.0, w_model(z_line_um, fit_x['w0'], fit_x['z0'], fit_x['M'], wavelength_um), label="fit")
        ax1.plot(z_um / 1000.0, widths_x, 'o', label="data")
        ax1.set_title("Profile waist in x")
        ax1.set_xlabel("z [mm]")
        ax1.set_ylabel("w(z) [µm]")
        ax1.legend()

        ax2.plot(z_line_um / 1000.0, w_model(z_line_um, fit_y['w0'], fit_y['z0'], fit_y['M'], wavelength_um), label="fit")
        ax2.plot(z_um / 1000.0, widths_y, 'o', label="data")
        ax2.set_title("Profile waist in y")
        ax2.set_xlabel("z [mm]")
        ax2.set_ylabel("w(z) [µm]")
        ax2.legend()

        layout.addWidget(canvas)
        toolbar = NavigationToolbar(canvas, dlg)
        layout.addWidget(toolbar)

        # Results table (wavelength, w0, z0, chi, M^2, Rayleigh length)
        tbl = QTableWidget(6, 2)
        tbl.setHorizontalHeaderLabels(["x waist", "y waist"])
        tbl.setVerticalHeaderLabels([
            "wavelength [nm]",
            "w0 [µm]",
            "z0 [mm]",
            "normalized chi-square",
            "M-square",
            "Rayleigh length [mm]"
        ])
        tbl.setItem(0, 0, QTableWidgetItem(f"{wavelength_um*1000:.6g}"))
        tbl.setItem(0, 1, QTableWidgetItem(f"{wavelength_um*1000:.6g}"))

        tbl.setItem(1, 0, QTableWidgetItem(f"{fit_x['w0']:.6g}"))
        tbl.setItem(1, 1, QTableWidgetItem(f"{fit_y['w0']:.6g}"))

        tbl.setItem(2, 0, QTableWidgetItem(f"{fit_x['z0']/1000.0:.6g}"))
        tbl.setItem(2, 1, QTableWidgetItem(f"{fit_y['z0']/1000.0:.6g}"))

        tbl.setItem(3, 0, QTableWidgetItem(f"{fit_x['chi']:.6g}"))
        tbl.setItem(3, 1, QTableWidgetItem(f"{fit_y['chi']:.6g}"))

        tbl.setItem(4, 0, QTableWidgetItem(f"{fit_x['M']:.6g}"))
        tbl.setItem(4, 1, QTableWidgetItem(f"{fit_y['M']:.6g}"))

        tbl.setItem(5, 0, QTableWidgetItem(f"{fit_x['zR']/1000.0:.6g}"))
        tbl.setItem(5, 1, QTableWidgetItem(f"{fit_y['zR']/1000.0:.6g}"))

        tbl.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(tbl)

        # Save button with full extended file (as in MATLAB pushbutton4_Callback)
        btn_save = QPushButton("Save Results (extended)")
        layout.addWidget(btn_save)

        def save_extended():
            if not self.entries:
                return
            # Use first image's basic values + profile outputs
            fr = self.entries[0].fit
            ps = self.pixel_size_um

            # Build line like MATLAB's extended file
            header = (
                "imagenumber normalized_chi-square waist_x waist_y I_max Offset Angle "
                "wavelength w_0_x w_0_y z_0_x z_0_y chi_square_px chi_square_py "
                "M_squared_x M_squared_y Rayleigh_length_x Rayleigh_length_y\n"
            )
            vals = [
                1, fr.chi, fr.wbx * ps, fr.wby * ps, fr.Imax, fr.offst, fr.angle_deg,
                wavelength_um * 1000.0,
                fit_x['w0'], fit_y['w0'],
                fit_x['z0'] / 1000.0, fit_y['z0'] / 1000.0,
                fit_x['chi'], fit_y['chi'],
                fit_x['M'], fit_y['M'],
                fit_x['zR'] / 1000.0, fit_y['zR'] / 1000.0
            ]

            path, _ = QFileDialog.getSaveFileName(self, "Save Extended Results", "beamFitFile1.txt", "Text (*.txt)")
            if not path:
                return
            with open(path, "w", encoding="utf-8") as f:
                f.write(header)
                f.write(
                    "{:6.0f} {:6.4f} {:6.4f} {:6.4f} {:6.4f} {:6.4f} {:6.4f} "
                    "{:6.4f} {:6.4f} {:6.4f} {:6.4f} {:6.4f} {:6.4f} {:6.4f} "
                    "{:6.4f} {:6.4f} {:6.4f} {:6.4f}\n".format(*vals)
                )
            QMessageBox.information(self, "Saved", f"Extended results saved to {path}")

        btn_save.clicked.connect(save_extended)

        dlg.resize(900, 900)
        dlg.exec_()

    # ------------------------------------------------------------------ #


def main():
    app = QApplication(sys.argv)
    win = GaussBeamApp()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()