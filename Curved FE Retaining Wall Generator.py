import sys
import os
import math
import time
import numpy as np

# --- STAAD.Pro connection (OpenSTAAD) ---
try:
    from openstaadpy import os_analytical
except ImportError:
    print("Warning: 'openstaadpy' not found. STAAD.Pro connection will fail.")
    os_analytical = None

# --- PySide6 GUI ---
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QTabWidget, QVBoxLayout, QHBoxLayout,
    QFormLayout, QLabel, QLineEdit, QPushButton, QGroupBox, QSplitter,
    QMessageBox, QComboBox, QCheckBox, QFileDialog
)
from PySide6.QtCore import Qt

# --- PyVista 3D viewer ---
import pyvista as pv
from pyvistaqt import QtInteractor


# ==================================================================
#  MATERIAL / SUPPORT CONSTANTS (SI: kN, m)
# ==================================================================
CONCRETE_MATERIAL_NAME = "CONCRETE"
CONCRETE_E_KN_M2 = 2.17185e7
CONCRETE_POISSON = 0.17
CONCRETE_DENSITY_KN_M3 = 23.5615
CONCRETE_ALPHA = 5.5e-6
CONCRETE_DAMPING = 0.05
CONCRETE_FCU_KN_M2 = 25000.0

DEFAULT_SOIL_SUBGRADE_KN_M3 = 10000.0
DEFAULT_LATERAL_SUBGRADE_KN_M3 = 5000.0

# OpenSTAAD support direction codes
DIRECTION_X = 0
DIRECTION_Y = 1
DIRECTION_Z = 2
SPRING_TYPE_COMPRESSION_ONLY = 1

# --- Inclined-support codes (used for curved-wall lateral soil springs) ---
INCLINED_TYPE_FIXEDBUT = 3
INCLINED_REFTYPE_DISTANCES = 0

# ------------------------------------------------------------------
# LOAD CASE NUMBERS (fixed layout)
#   LC1 = Dead load: Self Weight (+ soil surcharge on heel / toe, if entered)
#   LC2 = Water Hydrostatic
#   LC3 = Lateral Earth Pressure
#   LC4 = REPEAT: Dead (LC1) + Earth
#   LC5 = REPEAT: Dead (LC1) + Earth + Water
# ------------------------------------------------------------------
SELFWEIGHT_LOAD_CASE = 1
WATER_LOAD_CASE = 2
EARTH_LOAD_CASE = 3
REPEAT_LC_SW_EARTH = 4
REPEAT_LC_SW_EARTH_WATER = 5
LOAD_TYPE_DEAD = 0

# Soil surcharge is applied as a pressure in global Y (OpenSTAAD direction code 5 = GY),
# inside the dead load case (LC1), together with the self-weight.
SURCHARGE_DIRECTION_CODE = 5

# ------------------------------------------------------------------
# Curved-wall lateral soil method
# ------------------------------------------------------------------
CURVED_METHOD_MEMBERS = "Compression-Only Members"
CURVED_METHOD_SPRINGS = "Inclined Springs"
CURVED_LATERAL_METHOD_OPTIONS = [CURVED_METHOD_MEMBERS, CURVED_METHOD_SPRINGS]
DEFAULT_CURVED_LATERAL_METHOD = CURVED_METHOD_MEMBERS

# ------------------------------------------------------------------
# SOIL SIDE
# ------------------------------------------------------------------
SOIL_SIDE_AUTO = "Auto (toward curve center)"
SOIL_SIDE_RIGHT = "Right (+X)"
SOIL_SIDE_LEFT = "Left (-X)"
SOIL_SIDE_OPTIONS = [SOIL_SIDE_AUTO, SOIL_SIDE_RIGHT, SOIL_SIDE_LEFT]
DEFAULT_SOIL_SIDE = SOIL_SIDE_AUTO

# Default length of the equivalent soil strut (m).
DEFAULT_STRUT_LENGTH_M = 0.25

STRUT_REFERENCE_LENGTH_M = 1.0
STRUT_MATERIAL_NAME = CONCRETE_MATERIAL_NAME
STRUT_POISSON = CONCRETE_POISSON
STRUT_DENSITY_KN_M3 = CONCRETE_DENSITY_KN_M3
STRUT_E_KN_M2 = CONCRETE_E_KN_M2

# ------------------------------------------------------------------
# Compression-only handling for the CURVED-wall INCLINED-SPRING option only.
# ------------------------------------------------------------------
SPRING_TC_MODE = "COMPRESSION"
SPRING_TC_DIRECTION_TOKEN = "KFX"
SPRING_COMPRESSION_OUTPUT_FILENAME = "spring_compression_commands.txt"

# OpenSTAAD pressure direction codes
HYDRO_DIRECTION_OPTIONS = {
    "Local Z (Normal to Wall)": {"code": 3, "sign": -1},
    "Global X":                 {"code": 4, "sign": 1},
    "Global Y":                 {"code": 5, "sign": 1},
    "Global Z":                 {"code": 6, "sign": 1},
}
DEFAULT_HYDRO_DIRECTION = "Local Z (Normal to Wall)"

# ------------------------------------------------------------------
# Lateral soil support behaviour for the STEM
# ------------------------------------------------------------------
LATERAL_MODE_AUTO_NORMAL = "Auto (Normal to Wall)"
LATERAL_MODE_GLOBAL_X = "Global X (Plate Mat)"
LATERAL_SUPPORT_MODE_OPTIONS = [LATERAL_MODE_AUTO_NORMAL, LATERAL_MODE_GLOBAL_X]
DEFAULT_LATERAL_SUPPORT_MODE = LATERAL_MODE_AUTO_NORMAL

# ------------------------------------------------------------------
# Counterfort SIDE / placement constants
# ------------------------------------------------------------------
COUNTERFORT_SIDE_HEEL = "heel"
COUNTERFORT_SIDE_TOE = "toe"
COUNTERFORT_SIDE_ALTERNATE_HEEL_FIRST = "alternate (heel first)"
COUNTERFORT_SIDE_ALTERNATE_TOE_FIRST = "alternate (toe first)"
COUNTERFORT_SIDE_PATTERN_OPTIONS = [
    COUNTERFORT_SIDE_HEEL,
    COUNTERFORT_SIDE_TOE,
    COUNTERFORT_SIDE_ALTERNATE_HEEL_FIRST,
    COUNTERFORT_SIDE_ALTERNATE_TOE_FIRST,
]

COUNTERFORT_MODE_UNIFORM = "Uniform Spacing"
COUNTERFORT_MODE_MANUAL = "Manual Positions"

PERFORM_ANALYSIS_PRINT_OPTIONS = {
    "No Print":            0,
    "Print Load Data":     1,
    "Print Statics Check": 2,
    "Print Statics Load":  3,
    "Print Mode Shapes":   4,
    "Print Both":          5,
    "Print All":           6,
}
DEFAULT_PERFORM_ANALYSIS_PRINT_OPTION = "No Print"

# Earth-pressure defaults
DEFAULT_EARTH_KA = 0.333
DEFAULT_EARTH_GAMMA_KN_M3 = 18.0

# ------------------------------------------------------------------
# Save-to-folder settings
# NewSTAADFile unit codes: length 4 = Meter, force 5 = KiloNewton
# ------------------------------------------------------------------
DEFAULT_SAVE_FILENAME = "RetainingWall"
SAVE_LENGTH_UNIT_METER = 4
SAVE_FORCE_UNIT_KILONEWTON = 5
SAVE_WAIT_TIMEOUT_S = 20.0
INVALID_FILENAME_CHARS = '<>:"/\\|?*'

# Save methods + .std writer settings
SAVE_METHOD_DIRECT = "Fast: write .std directly (no live STAAD calls)"
SAVE_METHOD_LIVE = "Live: build inside STAAD.Pro, then SaveModel"
SAVE_METHOD_OPTIONS = [SAVE_METHOD_DIRECT, SAVE_METHOD_LIVE]
DEFAULT_SAVE_METHOD = SAVE_METHOD_DIRECT

STD_MAX_LINE = 79
STD_SUPPORT_DIR_NAMES = {0: "X", 1: "Y", 2: "Z"}
STD_AXIS_NAMES = {1: "X", 2: "Y", 3: "Z"}
STD_PRESSURE_DIR_NAMES = {3: "", 4: "GX", 5: "GY", 6: "GZ"}
STD_PRINT_OPTIONS = {
    0: "",
    1: " PRINT LOAD DATA",
    2: " PRINT STATICS CHECK",
    3: " PRINT STATICS LOAD",
    4: " PRINT MODE SHAPES",
    5: " PRINT BOTH",
    6: " PRINT ALL",
}


# ==================================================================
#  GEOMETRY / HELPER UTILITIES
# ==================================================================
def triangle_area(p1, p2, p3):
    return 0.5 * abs(p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))


def quad_area_3d(p1, p2, p3, p4):
    p1, p2, p3, p4 = (np.asarray(p, dtype=float) for p in (p1, p2, p3, p4))
    a1 = 0.5 * np.linalg.norm(np.cross(p2 - p1, p3 - p1))
    a2 = 0.5 * np.linalg.norm(np.cross(p3 - p1, p4 - p1))
    return float(a1 + a2)


def clean_polygon(poly_xy):
    if len(poly_xy) < 3:
        return []
    pts = []
    for pt in poly_xy:
        if not pts or (abs(pts[-1][0] - pt[0]) > 1e-4 or abs(pts[-1][1] - pt[1]) > 1e-4):
            pts.append(pt)
    if len(pts) > 1 and abs(pts[0][0] - pts[-1][0]) < 1e-4 and abs(pts[0][1] - pts[-1][1]) < 1e-4:
        pts.pop()
    if len(pts) < 3:
        return []
    cleaned = []
    n = len(pts)
    for i in range(n):
        prev_p = pts[(i - 1) % n]
        curr_p = pts[i]
        next_p = pts[(i + 1) % n]
        cross = (curr_p[0] - prev_p[0]) * (next_p[1] - prev_p[1]) - (curr_p[1] - prev_p[1]) * (next_p[0] - prev_p[0])
        if abs(cross) > 1e-5:
            cleaned.append(curr_p)
    return cleaned if len(cleaned) >= 3 else []


def compress_ids_to_ranges(ids):
    if not ids:
        return []
    ids = sorted(set(int(i) for i in ids))
    ranges = []
    start = prev = ids[0]
    rest = ids[1:]  # split out so the loop line has no bracket before the colon
    for x in rest:
        if x == prev + 1:
            prev = x
            continue
        ranges.append(f"{start} TO {prev}" if prev > start else f"{start}")
        start = prev = x
    ranges.append(f"{start} TO {prev}" if prev > start else f"{start}")
    return ranges


def build_spring_compression_block(node_ids, direction_token=SPRING_TC_DIRECTION_TOKEN,
                                   tc_mode=SPRING_TC_MODE, n_load_cases=2, per_line=8):
    ranges = compress_ids_to_ranges(node_ids)
    if not ranges:
        return ""

    lines = []
    for i in range(0, len(ranges), per_line):
        lines.append(" ".join(ranges[i:i + per_line]))

    spring_lines = []
    if len(lines) == 1:
        spring_lines.append(f"SPRING {tc_mode} {lines[0]} {direction_token}")
    else:
        spring_lines.append(f"SPRING {tc_mode} {lines[0]} -")
        middle = lines[1:-1]  # split out so the loop line has no bracket before the colon
        for mid in middle:
            spring_lines.append(f"    {mid} -")
        spring_lines.append(f"    {lines[-1]} {direction_token}")

    block = []
    block.append("*" + "=" * 66)
    block.append("* AUTO-GENERATED ONE-WAY (COMPRESSION-ONLY) SPRING COMMANDS")
    block.append("* OpenSTAAD cannot write these, so paste them into the STAAD editor.")
    block.append("*   1) Put 'SET NL {}' as the VERY FIRST line of the input file.".format(max(1, int(n_load_cases))))
    block.append("*   2) Put the 'SPRING {}' block below just BEFORE the load cases.".format(tc_mode))
    block.append("*   3) Ensure a 'CHANGE' follows the PERFORM ANALYSIS (already noted).")
    block.append("*   4) Optionally add 'SET ITERLIM 15' if convergence needs more cycles.")
    block.append("*" + "=" * 66)
    block.append(f"SET NL {max(1, int(n_load_cases))}")
    block.append("SET ITERLIM 15")
    block.append("")
    block.extend(spring_lines)
    block.append("")
    block.append("* ... your LOAD cases / PERFORM ANALYSIS here ...")
    block.append("CHANGE")
    block.append("*" + "=" * 66)
    return "\n".join(block)


# ------------------------------------------------------------------
#  SOIL-SIDE STRUT DIRECTION HELPER (shared by generator + preview)
# ------------------------------------------------------------------
def strut_direction(gx, gz, use_curve, cx, cz, soil_side):
    if use_curve:
        dx = cx - gx
        dz = cz - gz
        dist = math.hypot(dx, dz)
        if dist < 1e-12:
            nx, nz = 1.0, 0.0
        else:
            nx, nz = dx / dist, dz / dist
    else:
        dx = dz = 0.0
        dist = float("inf")
        nx, nz = 1.0, 0.0

    if soil_side == SOIL_SIDE_RIGHT:
        if nx < 0:
            nx, nz = -nx, -nz
    elif soil_side == SOIL_SIDE_LEFT:
        if nx > 0:
            nx, nz = -nx, -nz

    points_inward = use_curve and (nx * dx + nz * dz) > 0.0
    return nx, nz, dist, points_inward


# ------------------------------------------------------------------
#  COUNTERFORT PLACEMENT HELPERS
# ------------------------------------------------------------------
def compute_counterfort_positions_uniform(length_z, spacing, step_z):
    if not spacing or spacing <= 0:
        return []
    n_bays = max(1, round(length_z / spacing))
    if n_bays < 2:
        return []
    actual_spacing = length_z / n_bays
    positions, seen = [], set()
    for i in range(1, n_bays):
        snapped = round(round(i * actual_spacing / step_z) * step_z, 6)
        if 0 < snapped < length_z and snapped not in seen:
            seen.add(snapped)
            positions.append(snapped)
    return positions


def parse_manual_positions_with_side(positions_str, length_z, step_z=None):
    results = {}
    if not positions_str:
        return []
    for token in positions_str.split(','):
        token = token.strip()
        if not token:
            continue
        if ':' in token:
            z_part, side_part = token.split(':', 1)
        else:
            z_part, side_part = token, ""
        z_part = z_part.strip()
        side_part = side_part.strip().lower()
        try:
            z = float(z_part)
        except ValueError:
            continue
        if step_z and step_z > 0:
            z = round(round(z / step_z) * step_z, 6)
        if not (0 < z < length_z):
            continue
        side = None
        if side_part in ("heel", "h"):
            side = COUNTERFORT_SIDE_HEEL
        elif side_part in ("toe", "t"):
            side = COUNTERFORT_SIDE_TOE
        results[z] = side
    return sorted(results.items())


def resolve_counterfort_specs(mode, positions_str, side_pattern, length_z, spacing, step_z):
    is_alternating = side_pattern in (COUNTERFORT_SIDE_ALTERNATE_HEEL_FIRST, COUNTERFORT_SIDE_ALTERNATE_TOE_FIRST)
    alt_start_side = COUNTERFORT_SIDE_TOE if side_pattern == COUNTERFORT_SIDE_ALTERNATE_TOE_FIRST else COUNTERFORT_SIDE_HEEL

    def next_alt_side(idx):
        other = COUNTERFORT_SIDE_TOE if alt_start_side == COUNTERFORT_SIDE_HEEL else COUNTERFORT_SIDE_HEEL
        return alt_start_side if idx % 2 == 0 else other

    if mode == COUNTERFORT_MODE_MANUAL:
        raw = parse_manual_positions_with_side(positions_str, length_z, step_z)
        specs = []
        alt_idx = 0
        for z, side in raw:
            if side is None:
                if is_alternating:
                    side = next_alt_side(alt_idx)
                    alt_idx += 1
                else:
                    side = side_pattern
            specs.append((z, side))
        return specs
    else:
        zs = compute_counterfort_positions_uniform(length_z, spacing, step_z)
        specs = []
        for idx, z in enumerate(zs):
            side = next_alt_side(idx) if is_alternating else side_pattern
            specs.append((z, side))
        return specs


# ------------------------------------------------------------------
#  COUNTERFORT HEIGHT HELPER (shared by generator + preview)
#  requested <= 0 / None  -> full wall height (original behaviour)
#  otherwise              -> capped at wall height and snapped to the
#                            nearest stem mesh row so nodes stay connected
# ------------------------------------------------------------------
def snap_counterfort_height(requested, height_y, step_y=None):
    if not requested or requested <= 0:
        return float(height_y)
    ch = min(float(requested), float(height_y))
    if step_y and step_y > 0:
        ch = round(round(ch / step_y) * step_y, 6)
        ch = max(step_y, min(ch, float(height_y)))
    return float(ch)


# ------------------------------------------------------------------
#  SAVE-TO-FOLDER HELPERS
# ------------------------------------------------------------------
def build_model_path(folder, file_name):
    name = (file_name or "").strip()
    for ch in INVALID_FILENAME_CHARS:
        name = name.replace(ch, "_")
    if not name:
        name = DEFAULT_SAVE_FILENAME
    if not name.lower().endswith(".std"):
        name += ".std"
    return os.path.normpath(os.path.join(folder, name))


def make_unique_path(path):
    """Never overwrite an existing model: add _1, _2, ... if the name is taken."""
    if not os.path.exists(path):
        return path
    base, ext = os.path.splitext(path)
    i = 1
    while os.path.exists(f"{base}_{i}{ext}"):
        i += 1
    return f"{base}_{i}{ext}"


def _same_path(a, b):
    try:
        return os.path.normcase(os.path.abspath(str(a))) == os.path.normcase(os.path.abspath(str(b)))
    except Exception:
        return False


def _active_file_matches(staad, full_path):
    for args in ((True,), ()):
        try:
            cur = staad.GetSTAADFile(*args)
        except Exception:
            continue
        if not cur:
            continue
        cur = str(cur)
        if os.path.dirname(cur):
            return _same_path(cur, full_path)
        return cur.lower() == os.path.basename(full_path).lower()
    return False


def create_new_staad_file(staad, full_path, length_unit=SAVE_LENGTH_UNIT_METER,
                          force_unit=SAVE_FORCE_UNIT_KILONEWTON, timeout_s=SAVE_WAIT_TIMEOUT_S):
    """Create a new .STD file at full_path (Meter / kN) and make it the active model.
    Returns (staad, ok, error_message)."""
    try:
        staad.NewSTAADFile(full_path, int(length_unit), int(force_unit))
    except Exception as e:
        return staad, False, f"NewSTAADFile failed: {e}"

    # NewSTAADFile can return before the file is fully ready, so wait for it.
    t0 = time.time()
    ready = False
    while time.time() - t0 < timeout_s:
        if _active_file_matches(staad, full_path):
            ready = True
            break
        time.sleep(0.25)
    if not ready:
        return staad, False, "STAAD.Pro did not switch to the new file in time."
    time.sleep(0.5)

    try:
        fresh = os_analytical.connect()
        if fresh is not None:
            staad = fresh
    except Exception:
        pass
    return staad, True, None


def save_staad_model(staad):
    """Save the active model in place. Returns (ok, error_message)."""
    last_err = None
    for args in ((1,), ()):
        try:
            staad.SaveModel(*args)
            return True, None
        except Exception as e:
            last_err = e
    return False, str(last_err)


def open_std_in_staad(path):
    """Open a saved .std in the running STAAD.Pro instance (best effort)."""
    if not os_analytical:
        return False, "openstaadpy is not installed."
    try:
        staad = os_analytical.connect()
        staad.OpenSTAADFile(path)
        return True, None
    except Exception as e:
        return False, str(e)


# ==================================================================
#  IN-MEMORY MODEL RECORDER  ->  writes a STAAD .std file
#  Stands in for the OpenSTAAD root object: the generator below calls
#  the same Geometry / Property / Support / Load / Command methods, but
#  everything is stored in plain Python lists and dicts. One text file
#  is written at the end, so STAAD.Pro is never touched during the build.
# ==================================================================
def _std_num(v):
    """Coordinate text: 6 decimals, trailing zeros trimmed, never '-0'."""
    v = round(float(v), 6)
    if v == 0:
        return "0"
    return f"{v:.6f}".rstrip("0").rstrip(".")


def _std_g(v):
    """General number text for properties and loads."""
    return f"{float(v) + 0.0:.10g}"


def _as_list(x):
    if isinstance(x, (list, tuple, set)):
        return list(x)
    return [x]


def _wrap_tokens(tokens, width=STD_MAX_LINE):
    """Join tokens into lines no wider than the STAAD input width.
    Lines that continue end with ' -' (STAAD continuation)."""
    lines = []
    cur = ""
    for tok in tokens:
        if not cur:
            cur = tok
        elif len(cur) + 1 + len(tok) > width - 4:
            lines.append(cur + " -")
            cur = tok
        else:
            cur += " " + tok
    if cur:
        lines.append(cur)
    return lines


class _RecGeometry:
    def __init__(self, model):
        self.m = model

    def AddNode(self, x, y, z):
        self.m.nodes.append((float(x), float(y), float(z)))
        return len(self.m.nodes)

    def AddPlate(self, n1, n2, n3, n4=0):
        pid = self.m.next_elem()
        self.m.plates[pid] = (int(n1), int(n2), int(n3), int(n4) if n4 else 0)
        return pid

    def AddBeam(self, n1, n2):
        mid = self.m.next_elem()
        self.m.members[mid] = (int(n1), int(n2))
        return mid


class _RecProperty:
    def __init__(self, model):
        self.m = model

    def GetMaterialProperty(self, name):
        return self.m.materials.get(name)

    def CreateIsotropicMaterialConcrete(self, name, E, poisson, G, density, alpha, damp, fcu, flag=1):
        self.m.materials[name] = {"E": E, "poisson": poisson, "G": G, "density": density,
                                  "alpha": alpha, "damp": damp, "fcu": fcu}
        return True

    def CreateIsotropicMaterialProperties(self, name, E, poisson, G, density, alpha, damp):
        self.m.materials[name] = {"E": E, "poisson": poisson, "G": G, "density": density,
                                  "alpha": alpha, "damp": damp, "fcu": None}
        return True

    def CreatePlateThicknessProperty(self, thicknesses):
        pid = self.m.new_id()
        self.m.thick_props[pid] = float(_as_list(thicknesses)[0])
        return pid

    def AssignPlateThickness(self, plate_ids, prop_id):
        t = self.m.thick_props.get(prop_id)
        if t is None:
            return False
        for p in _as_list(plate_ids):
            self.m.plate_thick[int(p)] = t
        return True

    def AssignMaterialToPlate(self, name, plate_ids):
        for p in _as_list(plate_ids):
            self.m.plate_mat[int(p)] = name
        return True

    def CreatePrismaticCircleProperty(self, diameter):
        pid = self.m.new_id()
        self.m.dia_props[pid] = float(_as_list(diameter)[0])
        return pid

    def AssignBeamProperty(self, member_ids, prop_id):
        d = self.m.dia_props.get(prop_id)
        if d is None:
            return False
        for mid in _as_list(member_ids):
            self.m.member_dia[int(mid)] = d
        return True

    def AssignMaterialToMember(self, name, member_ids):
        for mid in _as_list(member_ids):
            self.m.member_mat[int(mid)] = name
        return True

    def GetBeamMaterialName(self, member_id):
        return self.m.member_mat.get(int(member_id), "")

    def CreateMemberCompressionSpec(self):
        pid = self.m.new_id()
        self.m.spec_props[pid] = "COMPRESSION"
        return pid

    def CreateMemberTrussSpec(self):
        pid = self.m.new_id()
        self.m.spec_props[pid] = "TRUSS"
        return pid

    def AssignMemberSpecToBeam(self, member_ids, spec_id):
        kind = self.m.spec_props.get(spec_id)
        if kind is None:
            return False
        for mid in _as_list(member_ids):
            self.m.member_spec[int(mid)] = kind
        return True


class _RecSupport:
    def __init__(self, model):
        self.m = model

    def CreatePlateMat(self, direction, subgrades, print_flag=0, spring_type=0):
        sid = self.m.new_id()
        k = float(_as_list(subgrades)[0])
        self.m.sup_defs[sid] = ("platemat", int(direction), k, int(spring_type))
        return sid

    def CreateSupportFixed(self):
        sid = self.m.new_id()
        self.m.sup_defs[sid] = ("fixed",)
        return sid

    def CreateInclinedSupport(self, inc_type, ref_type, ref_node, coord, release_spec, spring_spec):
        sid = self.m.new_id()
        self.m.sup_defs[sid] = ("inclined", int(inc_type), int(ref_type), int(ref_node),
                                tuple(coord), tuple(release_spec), tuple(spring_spec))
        return sid

    def AssignSupportToEntityList(self, support_id, ids):
        bucket = self.m.sup_assign.setdefault(support_id, set())
        for i in _as_list(ids):
            bucket.add(int(i))
        return True

    def AssignSupportToNode(self, node, support_id):
        self.m.sup_assign.setdefault(support_id, set()).add(int(node))
        return True


class _RecLoad:
    def __init__(self, model):
        self.m = model

    def SetLoadActive(self, lc):
        lc = int(lc)
        if lc in self.m.loads:
            self.m.active_load = lc
            return True
        return False

    def CreateNewPrimaryLoadEx2(self, title, load_type, lc):
        lc = int(lc)
        if lc not in self.m.loads:
            self.m.loads[lc] = {"title": str(title), "selfweight": [], "pressure": [], "repeat": [],
                                "trap": []}  # linearly varying plate pressures
        self.m.active_load = lc
        return True

    def _active(self):
        return self.m.loads[self.m.active_load]

    def AddSelfWeightInXYZ(self, direction, factor):
        self._active()["selfweight"].append((int(direction), float(factor)))
        return True

    def AddElementHydrostaticPressure(self, *args, **kwargs):
        # Not written natively: the generator then falls back to per-plate pressures.
        raise NotImplementedError("per-plate pressures are used when writing a .std file")

    def AddElementPressure(self, plate_ids, direction, value, *extra):
        for p in _as_list(plate_ids):
            self._active()["pressure"].append((int(p), int(direction), float(value)))
        return True

    def AddRepeatLoad(self, lc_list, factor_list):
        self._active()["repeat"] = [(int(a), float(f)) for a, f in zip(lc_list, factor_list)]
        return True

    # Linearly varying pressure on a plate: value_bottom at the plate's lower edge,
    # value_top at its upper edge. Written to the .std as
    #   element-list TRAP {GX|GY|GZ} Y f1 f2     (direction omitted = normal to plate, local Z)
    # which STAAD.Pro's plate "Trapezoidal" load uses. Local Y runs up the wall for the stem plates
    # (local x = node1->node2 = along the wall, z = IJ x JL, y = z x x).
    def AddPlateLinearPressure(self, plate_ids, direction, value_bottom, value_top):
        for p in _as_list(plate_ids):
            self._active()["trap"].append((int(p), int(direction), float(value_bottom), float(value_top)))
        return True


class _RecCommand:
    def __init__(self, model):
        self.m = model

    def PerformAnalysis(self, print_option=0):
        self.m.analysis_done = True
        self.m.analysis_print = int(print_option)
        return True


class StdModelRecorder:
    """In-memory model with an OpenSTAAD-like API that writes a .std file."""

    def __init__(self, title="Retaining Wall"):
        self.title = title
        self.nodes = []
        self.plates = {}
        self.members = {}
        self._elem_counter = 0
        self._id_counter = 0
        self.materials = {}
        self.thick_props = {}
        self.dia_props = {}
        self.spec_props = {}
        self.plate_thick = {}
        self.plate_mat = {}
        self.member_dia = {}
        self.member_mat = {}
        self.member_spec = {}
        self.sup_defs = {}
        self.sup_assign = {}
        self.loads = {}
        self.active_load = None
        self.analysis_done = False
        self.analysis_print = 0
        self.spring_compression_nodes = []
        self.Geometry = _RecGeometry(self)
        self.Property = _RecProperty(self)
        self.Support = _RecSupport(self)
        self.Load = _RecLoad(self)
        self.Command = _RecCommand(self)

    def next_elem(self):
        # Members and plates share one numbering so their numbers never clash.
        self._elem_counter += 1
        return self._elem_counter

    def new_id(self):
        self._id_counter += 1
        return self._id_counter

    def AddSpringCompression(self, node_ids):
        self.spring_compression_nodes = sorted(set(int(n) for n in node_ids))

    def summary(self):
        return f"{len(self.nodes)} nodes, {len(self.plates)} plates, {len(self.members)} members"

    # ---------------- text builders ----------------
    def _geometry_lines(self):
        out = ["JOINT COORDINATES"]
        for i, (x, y, z) in enumerate(self.nodes, start=1):
            out.append(f"{i} {_std_num(x)} {_std_num(y)} {_std_num(z)};")
        if self.members:
            out.append("MEMBER INCIDENCES")
            for mid in sorted(self.members):
                n1, n2 = self.members[mid]
                out.append(f"{mid} {n1} {n2};")
        if self.plates:
            out.append("ELEMENT INCIDENCES SHELL")
            for pid in sorted(self.plates):
                n1, n2, n3, n4 = self.plates[pid]
                tail = f"{n1} {n2} {n3} {n4}" if n4 else f"{n1} {n2} {n3}"
                out.append(f"{pid} {tail};")
        return out

    def _property_lines(self):
        out = []
        if self.materials:
            out.append("DEFINE MATERIAL START")
            for name, p in self.materials.items():
                out.append(f"ISOTROPIC {name}")
                out.append(f"E {_std_g(p['E'])}")
                out.append(f"POISSON {_std_g(p['poisson'])}")
                out.append(f"DENSITY {_std_g(p['density'])}")
                out.append(f"ALPHA {_std_g(p['alpha'])}")
                out.append(f"DAMP {_std_g(p['damp'])}")
                out.append(f"G {_std_g(p['G'])}")
                if p.get("fcu"):
                    out.append("TYPE CONCRETE")
                    out.append(f"STRENGTH FCU {_std_g(p['fcu'])}")
            out.append("END DEFINE MATERIAL")

        by_thick = {}
        for pid, t in self.plate_thick.items():
            by_thick.setdefault(round(t, 6), []).append(pid)
        if by_thick:
            out.append("ELEMENT PROPERTY")
            for t in sorted(by_thick):
                out.extend(_wrap_tokens(compress_ids_to_ranges(by_thick[t]) + ["THICKNESS " + _std_g(t)]))

        by_dia = {}
        for mid, d in self.member_dia.items():
            by_dia.setdefault(round(d, 8), []).append(mid)
        if by_dia:
            out.append("MEMBER PROPERTY")
            for d in sorted(by_dia):
                out.extend(_wrap_tokens(compress_ids_to_ranges(by_dia[d]) + ["PRIS YD " + _std_g(d)]))

        names = set(self.plate_mat.values()) | set(self.member_mat.values())
        if names:
            out.append("CONSTANTS")
            for name in sorted(names):  # this generator only ever uses one material
                out.append(f"MATERIAL {name} ALL")

        for kind in ("TRUSS", "COMPRESSION"):
            ids = sorted(m for m, k in self.member_spec.items() if k == kind)
            if ids:
                out.append("MEMBER " + kind)
                out.extend(_wrap_tokens(compress_ids_to_ranges(ids)))
        return out

    def _support_lines(self):
        body = []
        for sid, nodes in self.sup_assign.items():
            d = self.sup_defs.get(sid)
            if not d or not nodes:
                continue
            ids = compress_ids_to_ranges(sorted(nodes))
            kind = d[0]
            if kind == "platemat":
                tail = "PLATE MAT DIRECTION " + STD_SUPPORT_DIR_NAMES.get(d[1], "Y") + " SUBGRADE " + _std_g(d[2])
                if d[3] == SPRING_TYPE_COMPRESSION_ONLY:
                    tail += " COMP"
                body.extend(_wrap_tokens(ids + [tail]))
            elif kind == "fixed":
                body.extend(_wrap_tokens(ids + ["FIXED"]))
            elif kind == "inclined":
                inc_type, ref_type, ref_node, coord, rel, spr = d[1], d[2], d[3], d[4], d[5], d[6]
                if ref_type == 1:
                    ref = "REF " + " ".join(_std_g(c) for c in coord)
                elif ref_type == 2:
                    ref = "REFJT " + str(ref_node)
                else:
                    ref = " ".join(_std_g(c) for c in coord)
                words = ["INCLINED"] + ref.split()
                if inc_type == 1:
                    words.append("PINNED")
                elif inc_type == 2:
                    words.append("FIXED")
                else:
                    dof = ["FX", "FY", "FZ", "MX", "MY", "MZ"]
                    released = [dof[i] for i in range(6) if rel[i] == 1]
                    springs = []
                    for i in range(6):
                        if rel[i] == -1:
                            springs.extend(["K" + dof[i], _std_g(spr[i])])
                    words.extend(["FIXED", "BUT"] + released + springs)
                body.extend(_wrap_tokens(ids + words))
        out = []
        if body:
            out.append("SUPPORTS")
            out.extend(body)
        if self.spring_compression_nodes:
            ids = compress_ids_to_ranges(self.spring_compression_nodes)
            out.extend(_wrap_tokens(["SPRING COMPRESSION"] + ids + [SPRING_TC_DIRECTION_TOKEN]))
        return out

    def _load_lines(self):
        out = []
        for lc, ld in self.loads.items():
            out.append(f"LOAD {lc} TITLE {ld['title']}")
            repeat_pairs = ld["repeat"]
            selfweights = ld["selfweight"]
            pressures = ld["pressure"]
            traps = ld["trap"]
            if repeat_pairs:
                out.append("REPEAT LOAD")
                out.extend(_wrap_tokens([f"{a} {_std_g(f)}" for a, f in repeat_pairs]))
            for direction, factor in selfweights:
                out.append(f"SELFWEIGHT {STD_AXIS_NAMES.get(direction, 'Y')} {_std_g(factor)}")
            if pressures or traps:
                out.append("ELEMENT LOAD")
            if pressures:
                groups = {}
                for pid, direction, val in pressures:
                    groups.setdefault((direction, round(val, 6)), []).append(pid)
                for key in sorted(groups):
                    dname = STD_PRESSURE_DIR_NAMES.get(key[0], "")
                    tail = "PRESSURE " + (dname + " " if dname else "") + _std_g(key[1])
                    out.extend(_wrap_tokens(compress_ids_to_ranges(groups[key]) + [tail]))
            if traps:
                tgroups = {}
                for pid, direction, f_bot, f_top in traps:
                    tgroups.setdefault((direction, round(f_bot, 6), round(f_top, 6)), []).append(pid)
                for key in sorted(tgroups):
                    dname = STD_PRESSURE_DIR_NAMES.get(key[0], "")
                    tail = ("TRAP " + (dname + " " if dname else "") + "Y "
                            + _std_g(key[1]) + " " + _std_g(key[2]))
                    out.extend(_wrap_tokens(compress_ids_to_ranges(tgroups[key]) + [tail]))
        return out

    def build_lines(self):
        lines = ["STAAD SPACE " + self.title,
                 "START JOB INFORMATION",
                 "ENGINEER DATE " + time.strftime("%d-%b-%y"),
                 "END JOB INFORMATION"]
        if self.spring_compression_nodes:
            lines.append("SET NL " + str(max(1, len(self.loads))))
            lines.append("SET ITERLIM 15")
        lines.append("INPUT WIDTH " + str(STD_MAX_LINE))
        lines.append("UNIT METER KN")
        lines.extend(self._geometry_lines())
        lines.extend(self._property_lines())
        lines.extend(self._support_lines())
        lines.extend(self._load_lines())
        if self.analysis_done:
            lines.append("PERFORM ANALYSIS" + STD_PRINT_OPTIONS.get(self.analysis_print, ""))
            if self.spring_compression_nodes:
                lines.append("CHANGE")
        lines.append("FINISH")
        return lines

    def write(self, path):
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(self.build_lines()) + "\n")
        return path


# ==================================================================
#  OPENSTAAD BACKEND FUNCTIONS
# ==================================================================
def ensure_concrete_material(prop, material_name=CONCRETE_MATERIAL_NAME):
    try:
        existing = prop.GetMaterialProperty(material_name)
        if existing:
            return True
    except Exception:
        pass
    try:
        shear_modulus = CONCRETE_E_KN_M2 / (2.0 * (1.0 + CONCRETE_POISSON))
        prop.CreateIsotropicMaterialConcrete(
            material_name, CONCRETE_E_KN_M2, CONCRETE_POISSON, shear_modulus,
            CONCRETE_DENSITY_KN_M3, CONCRETE_ALPHA, CONCRETE_DAMPING, CONCRETE_FCU_KN_M2, 1
        )
        return True
    except Exception as e:
        print(f"Error creating material '{material_name}': {e}")
        return False


def _material_exists(prop, name):
    try:
        res = prop.GetMaterialProperty(name)
        if res in (None, False):
            return False
        if isinstance(res, (list, tuple, str)) and len(res) == 0:
            return False
        return True
    except Exception:
        return False


def resolve_strut_material(prop, desired_name=STRUT_MATERIAL_NAME, desired_E=STRUT_E_KN_M2,
                           poisson=STRUT_POISSON, density=STRUT_DENSITY_KN_M3):
    if _material_exists(prop, desired_name):
        print(f"Strut material '{desired_name}' already defined; using it (E={desired_E:.4g}).")
        return desired_name, float(desired_E)

    G = desired_E / (2.0 * (1.0 + poisson))
    dens = float(density) if density and density > 0 else 1.0e-9

    create_attempts = [
        ("CreateIsotropicMaterialProperties", (desired_name, float(desired_E), float(poisson), float(G),
                                               dens, 0.0, 0.0)),
        ("CreateIsotropicMaterialConcrete",   (desired_name, float(desired_E), float(poisson), float(G),
                                               dens, 0.0, 0.0, float(CONCRETE_FCU_KN_M2), 1)),
    ]
    for method_name, args in create_attempts:
        method = getattr(prop, method_name, None)
        if method is None:
            print(f"  {method_name} not available on this build.")
            continue
        try:
            method(*args)
            print(f"Strut material '{desired_name}' created via {method_name} "
                  f"(E={desired_E:.4g}, density={dens:g}).")
            return desired_name, float(desired_E)
        except Exception as e:
            print(f"  {method_name}('{desired_name}') failed: {e}")
            continue

    ensure_concrete_material(prop, CONCRETE_MATERIAL_NAME)
    print(f"Dedicated strut material could not be created; falling back to "
          f"'{CONCRETE_MATERIAL_NAME}' for the struts (E={CONCRETE_E_KN_M2:.4g}).")
    return CONCRETE_MATERIAL_NAME, float(CONCRETE_E_KN_M2)


def _assign_plate_mat_support(sup, support_id, plate_ids, label="plate support"):
    plate_ids = [int(p) for p in plate_ids]
    if not plate_ids:
        return 0
    try:
        sup.AssignSupportToEntityList(support_id, plate_ids)
    except Exception:
        pass
    ok = 0
    for pid in plate_ids:
        try:
            sup.AssignSupportToEntityList(support_id, [pid])
            ok += 1
        except Exception:
            pass
    print(f"  {label}: assigned to {ok}/{len(plate_ids)} plates.")
    return ok


def fix_nodes(sup, node_ids, label="fixed support"):
    """Apply a FIXED support to NODES using the OFFICIAL signature
    AssignSupportToNode(node, support_id)."""
    node_ids = [int(n) for n in node_ids]
    if not node_ids:
        return 0

    try:
        fixed_id = sup.CreateSupportFixed()
    except Exception as e:
        print(f"  {label}: CreateSupportFixed failed: {e}")
        return 0

    node_assign = getattr(sup, "AssignSupportToNode", None)
    nodes_assign = getattr(sup, "AssignSupportToNodes", None)

    ok = 0
    for nid in node_ids:
        done = False
        if node_assign is not None:
            try:
                node_assign(nid, fixed_id)
                done = True
            except Exception:
                if nodes_assign is not None:
                    try:
                        nodes_assign([nid], fixed_id)
                        done = True
                    except Exception:
                        done = False
        if not done and node_assign is not None:
            try:
                node_assign(fixed_id, nid)
                done = True
            except Exception:
                done = False
        if not done:
            try:
                sup.AssignSupportToEntityList(fixed_id, [nid])
                done = True
            except Exception:
                done = False
        if done:
            ok += 1
        else:
            print(f"  {label}: could not fix node {nid}.")

    print(f"  {label}: assigned to {ok}/{len(node_ids)} nodes.")
    return ok


def add_soil_support(staad, base_plate_ids, subgrade_modulus):
    if not base_plate_ids:
        return
    try:
        sup = staad.Support
        try:
            plate_mat_id = sup.CreatePlateMat(DIRECTION_Y, [subgrade_modulus], 0, SPRING_TYPE_COMPRESSION_ONLY)
        except Exception:
            plate_mat_id = sup.CreatePlateMat(DIRECTION_Y, subgrade_modulus, 0, SPRING_TYPE_COMPRESSION_ONLY)
        _assign_plate_mat_support(sup, plate_mat_id, base_plate_ids, label="Base soil (plate mat)")
    except Exception as e:
        print(f"Support error: {e}")


def add_lateral_soil_support(staad, stem_plate_ids, subgrade_modulus, direction=DIRECTION_X,
                             compression_only=True):
    if not stem_plate_ids:
        return
    spring_type = SPRING_TYPE_COMPRESSION_ONLY if compression_only else 0
    try:
        sup = staad.Support
        try:
            plate_mat_id = sup.CreatePlateMat(direction, [subgrade_modulus], 0, spring_type)
        except Exception:
            plate_mat_id = sup.CreatePlateMat(direction, subgrade_modulus, 0, spring_type)
        _assign_plate_mat_support(sup, plate_mat_id, stem_plate_ids, label="Lateral soil (plate mat)")
    except Exception as e:
        print(f"Lateral Support error: {e}")


def assign_material_to_members(prop, material_name, member_ids):
    if not member_ids:
        return 0

    member_ids = [int(m) for m in member_ids]

    assigned = 0
    for mid in member_ids:
        ok = False
        for arg in (mid, [mid]):
            try:
                prop.AssignMaterialToMember(material_name, arg)
                ok = True
                break
            except Exception:
                continue
        if ok:
            assigned += 1
        else:
            print(f"  Could not assign material '{material_name}' to strut {mid}.")

    getter = getattr(prop, "GetBeamMaterialName", None)
    if getter is not None:
        def _confirmed_ids():
            good = []
            for mid in member_ids:
                try:
                    nm = getter(mid)
                    if nm and material_name.upper() in str(nm).upper():
                        good.append(mid)
                except Exception:
                    return None
            return good

        good = _confirmed_ids()
        if good is not None:
            missing = [m for m in member_ids if m not in set(good)]
            for mid in missing:
                for arg in (mid, [mid]):
                    try:
                        prop.AssignMaterialToMember(material_name, arg)
                        break
                    except Exception:
                        continue
            good2 = _confirmed_ids()
            confirmed = len(good2) if good2 is not None else assigned
            print(f"Material '{material_name}' confirmed on {confirmed}/{len(member_ids)} struts.")
            return confirmed

    print(f"Material '{material_name}' assigned to {assigned}/{len(member_ids)} struts (unverified).")
    return assigned


def add_lateral_soil_support_curved_members(staad, stem_node_trib, stem_node_coords,
                                            curve_radius, curve_turn_dir, stem_x,
                                            subgrade_modulus, strut_length=DEFAULT_STRUT_LENGTH_M,
                                            strut_E=STRUT_E_KN_M2, compression_only=True,
                                            strut_material=STRUT_MATERIAL_NAME,
                                            soil_side=DEFAULT_SOIL_SIDE):
    if not stem_node_trib:
        return 0, 0

    geom = staad.Geometry
    prop = staad.Property
    sup = staad.Support

    use_curve = bool(curve_radius) and curve_radius > 0
    cx = stem_x + curve_turn_dir * curve_radius if use_curve else 0.0
    cz = 0.0

    eff_material, eff_E = resolve_strut_material(prop, strut_material, strut_E)
    strut_E = float(eff_E)
    strut_material = eff_material

    axial_spec_id = None
    axial_spec_label = None
    if compression_only:
        try:
            axial_spec_id = prop.CreateMemberCompressionSpec()
            axial_spec_label = "MEMBER COMPRESSION"
        except Exception as e:
            print(f"CreateMemberCompressionSpec failed: {e}")
    if axial_spec_id is None:
        try:
            axial_spec_id = prop.CreateMemberTrussSpec()
            axial_spec_label = "MEMBER TRUSS"
            if compression_only:
                print("  -> MEMBER COMPRESSION unavailable; using MEMBER TRUSS instead.")
        except Exception as e:
            print(f"CreateMemberTrussSpec failed: {e}")

    member_ids = []
    anchor_node_ids = []
    area_buckets = {}

    skipped = 0
    too_long = 0
    for nid, area in stem_node_trib.items():
        if area <= 0:
            skipped += 1
            continue
        gx, gy, gz = stem_node_coords[nid]

        nx, nz, dist_axis, points_inward = strut_direction(gx, gz, use_curve, cx, cz, soil_side)

        L = float(strut_length)
        if points_inward and L >= dist_axis:
            L = 0.5 * dist_axis
            too_long += 1

        ax = gx + nx * L
        az = gz + nz * L
        try:
            anchor_id = int(geom.AddNode(float(ax), float(gy), float(az)))
        except Exception as e:
            print(f"AddNode (anchor) failed for stem node {nid}: {e}")
            continue
        try:
            mid = int(geom.AddBeam(int(nid), anchor_id))
        except Exception as e:
            print(f"AddBeam (strut) failed for stem node {nid}: {e}")
            continue

        K = float(subgrade_modulus) * float(area)
        A = K * L / float(strut_E)
        key = round(A, 10)
        area_buckets.setdefault(key, {"A": A, "members": []})["members"].append(mid)

        member_ids.append(mid)
        anchor_node_ids.append(anchor_id)

    if not member_ids:
        print("Lateral soil (members): no struts created.")
        return 0, 0

    for key, b in area_buckets.items():
        A = b["A"]
        d = math.sqrt(4.0 * A / math.pi)
        prop_id = None
        try:
            prop_id = prop.CreatePrismaticCircleProperty(float(d))
        except Exception:
            try:
                prop_id = prop.CreatePrismaticCircleProperty([float(d)])
            except Exception as e:
                print(f"CreatePrismaticCircleProperty failed (d={d:.5f}): {e}")
        if prop_id is None:
            continue
        sect_ok = 0
        members_here = b["members"]  # split out so the loop line has no bracket before the colon
        for mid in members_here:
            done = False
            for arg in (int(mid), [int(mid)]):
                try:
                    prop.AssignBeamProperty(arg, prop_id)
                    done = True
                    break
                except Exception:
                    continue
            if done:
                sect_ok += 1
            else:
                print(f"AssignBeamProperty failed for strut {mid}.")
        print(f"Section (d={d*1000:.2f} mm) assigned to {sect_ok}/{len(members_here)} struts.")

    assign_material_to_members(prop, strut_material, member_ids)

    def _assign_spec(spec_id, label):
        if spec_id is None:
            return
        ok = 0
        for mid in member_ids:
            done = False
            for arg in (int(mid), [int(mid)]):
                try:
                    prop.AssignMemberSpecToBeam(arg, spec_id)
                    done = True
                    break
                except Exception:
                    continue
            if done:
                ok += 1
            else:
                print(f"Assign {label} spec failed for strut {mid}.")
        print(f"{label} spec assigned to {ok}/{len(member_ids)} struts.")

    _assign_spec(axial_spec_id, axial_spec_label or "AXIAL")

    n_fixed = fix_nodes(sup, anchor_node_ids, label="Strut anchors (fixed)")
    if n_fixed < len(anchor_node_ids):
        print(f"  WARNING: only {n_fixed}/{len(anchor_node_ids)} strut anchors fixed - "
              "struts without a fixed anchor give NO lateral restraint.")

    note = f" [{axial_spec_label}]" if axial_spec_label else ""
    warn = f"  [{too_long} strut(s) shortened to fit inside curve axis]" if too_long else ""
    print(f"Lateral soil (members, soil_side={soil_side}): {len(member_ids)} struts{note} "
          f"created across {len(area_buckets)} section size(s).{warn}")
    return len(member_ids), len(member_ids)


def add_lateral_soil_support_curved_normal(staad, stem_node_trib, stem_node_coords,
                                           curve_radius, curve_turn_dir, stem_x,
                                           subgrade_modulus):
    if not stem_node_trib:
        return 0, []

    sup = staad.Support
    cx = stem_x + curve_turn_dir * curve_radius
    cz = 0.0

    buckets = {}
    skipped = 0
    for nid, area in stem_node_trib.items():
        if area <= 0:
            skipped += 1
            continue
        gx, gy, gz = stem_node_coords[nid]
        dx = cx - gx
        dz = cz - gz
        if abs(dx) < 1e-9 and abs(dz) < 1e-9:
            skipped += 1
            continue
        k = float(subgrade_modulus) * float(area)
        key = (round(dx, 5), round(dz, 5), round(k, 4))
        b = buckets.get(key)
        if b is None:
            buckets[key] = {"dir": (dx, dz), "k": k, "nodes": [int(nid)]}
        else:
            b["nodes"].append(int(nid))

    total_nodes = 0
    total_defs = 0
    unique_axes = set()
    sprung_node_ids = []
    node_assign = getattr(sup, "AssignSupportToNode", None)
    for key, b in buckets.items():
        dx, dz = b["dir"]
        k = b["k"]
        unique_axes.add((round(dx, 5), round(dz, 5)))
        release_spec = [-1, 1, 1, 1, 1, 1]
        spring_spec = [k, 0.0, 0.0, 0.0, 0.0, 0.0]
        try:
            sup_id = sup.CreateInclinedSupport(
                INCLINED_TYPE_FIXEDBUT, INCLINED_REFTYPE_DISTANCES, 0,
                [float(dx), 0.0, float(dz)], release_spec, spring_spec
            )
        except Exception as e:
            print(f"Inclined support create failed (axis={b['dir']}, k={k:.1f}): {e}")
            continue

        node_list = b["nodes"]
        assigned_here = 0
        for nid in node_list:
            done = False
            if node_assign is not None:
                try:
                    node_assign(int(nid), sup_id)
                    done = True
                except Exception:
                    try:
                        node_assign(sup_id, int(nid))
                        done = True
                    except Exception:
                        done = False
            if not done:
                try:
                    sup.AssignSupportToEntityList(sup_id, [int(nid)])
                    done = True
                except Exception:
                    done = False
            if done:
                assigned_here += 1

        if assigned_here > 0:
            total_defs += 1
            total_nodes += assigned_here
            sprung_node_ids.extend(node_list)

    sprung_node_ids = sorted(set(sprung_node_ids))
    print(f"Curved-wall lateral soil: {total_nodes} nodes sprung across "
          f"{total_defs} support definitions ({len(unique_axes)} distinct radial axes; "
          f"{skipped} nodes skipped).")
    return total_nodes, sprung_node_ids


def apply_or_emit_spring_compression(staad, node_ids, n_load_cases=2,
                                     tc_mode=SPRING_TC_MODE,
                                     direction_token=SPRING_TC_DIRECTION_TOKEN,
                                     output_dir=None):
    if not node_ids:
        return False, "", None

    applied_via_api = False

    command_block = build_spring_compression_block(
        node_ids, direction_token=direction_token, tc_mode=tc_mode, n_load_cases=n_load_cases
    )

    file_path = None
    try:
        target_dir = output_dir or os.getcwd()
        file_path = os.path.join(target_dir, SPRING_COMPRESSION_OUTPUT_FILENAME)
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(command_block + "\n")
        print(f"SPRING {tc_mode} command block written to: {file_path}")
    except Exception as e:
        print(f"Could not write spring-compression command file: {e}")
        file_path = None

    return applied_via_api, command_block, file_path


def add_selfweight_load(staad, load_case_number, title="Self Weight"):
    try:
        load_interface = staad.Load
        if not load_interface.SetLoadActive(load_case_number):
            load_interface.CreateNewPrimaryLoadEx2(str(title), 1, load_case_number)
            load_interface.SetLoadActive(load_case_number)
        load_interface.AddSelfWeightInXYZ(2, -1.0)
    except Exception as e:
        print(f"Self-Weight error: {e}")


def add_repeat_load_case(staad, lc_num, title, load_factor_pairs, load_type=LOAD_TYPE_DEAD):
    """
    Create a NEW PRIMARY load case (number lc_num, given title) that REPEATs prior
    primary load cases with factors - the documented OpenSTAAD pattern:
        Load.CreateNewPrimaryLoadEx2(title, type, lcNum)
        Load.AddRepeatLoad(loadCaseList, factorList)
    A REPEAT LOAD is a true primary case (correct for P-Delta), not a linear
    combination.

    load_factor_pairs : list of (primary_lc_no, factor). Pairs with factor None
    or 0 are dropped. If nothing remains, the case is not created.

    Returns True on success.
    """
    pairs = [(int(lc), float(f)) for (lc, f) in load_factor_pairs
             if f is not None and abs(float(f)) > 1e-12]
    if not pairs:
        print(f"Repeat load '{title}' (LC{lc_num}): no valid constituents; skipped.")
        return False

    try:
        load = staad.Load
        created = False
        try:
            load.CreateNewPrimaryLoadEx2(str(title), int(load_type), int(lc_num))
            created = True
        except Exception:
            try:
                load.CreateNewPrimaryLoadEx(str(title), int(load_type))
                created = True
            except Exception:
                try:
                    load.CreateNewPrimaryLoad(str(title))
                    created = True
                except Exception as e:
                    print(f"Could not create repeat load case '{title}': {e}")
        try:
            load.SetLoadActive(int(lc_num))
        except Exception:
            pass

        lc_list = [p[0] for p in pairs]
        f_list = [p[1] for p in pairs]

        ok = False
        for args in ((lc_list, f_list), (tuple(lc_list), tuple(f_list))):
            try:
                load.AddRepeatLoad(*args)
                ok = True
                break
            except Exception:
                continue
        if not ok:
            print(f"AddRepeatLoad failed for '{title}' (LC{lc_num}).")

        desc = " + ".join(f"{f:g}xLC{lc}" for lc, f in pairs)
        print(f"Repeat load '{title}' (LC{lc_num}) = {desc}  [{'OK' if ok and created else 'CHECK'}]")
        return ok
    except Exception as e:
        print(f"Repeat load error '{title}': {e}")
        return False


def apply_triangular_pressure_load(staad, stem_plates_for_load, load_case, load_at_top, load_at_bottom,
                                   pressure_height, direction_code=3, direction_sign=-1,
                                   load_name="Hydrostatic Load"):
    """Triangular (depth-increasing) pressure over the stem, used for BOTH water
    hydrostatic (gamma_water) and lateral earth pressure (Ka*gamma_soil)."""
    if not stem_plates_for_load:
        return
    try:
        load_interface = staad.Load
        try:
            load_interface.CreateNewPrimaryLoadEx2(str(load_name), 1, int(load_case))
        except Exception:
            pass
        load_interface.SetLoadActive(int(load_case))

        loaded_plate_ids = [int(pid) for (pid, y1, y2) in stem_plates_for_load if y2 <= pressure_height + 1e-6]
        if not loaded_plate_ids:
            return

        print(f"Applying '{load_name}' to {len(loaded_plate_ids)} plates in Load Case {load_case}...")

        # .std recorder only: write a smooth linear (TRAP) pressure on every plate instead of
        # one uniform value per plate, so the load follows the true hydrostatic gradient.
        gradient_fn = getattr(load_interface, "AddPlateLinearPressure", None)
        if gradient_fn is not None:
            def _p_at(y):
                yy = min(max(y, 0.0), pressure_height)
                return load_at_top + (load_at_bottom - load_at_top) * ((pressure_height - yy) / pressure_height)

            n_linear = 0
            for pid, y1, y2 in stem_plates_for_load:
                if y2 > pressure_height + 1e-6:
                    continue
                gradient_fn([int(pid)], direction_code,
                            float(direction_sign * _p_at(y1)), float(direction_sign * _p_at(y2)))
                n_linear += 1
            print(f"'{load_name}': linear (TRAP) pressure written on {n_linear} plates.")
            return

        applied_natively = False
        try:
            success = load_interface.AddElementHydrostaticPressure(
                loaded_plate_ids, direction_code, 2,
                float(direction_sign * load_at_bottom), float(direction_sign * load_at_top)
            )
            if success:
                applied_natively = True
        except Exception:
            try:
                success = load_interface.AddElementHydrostaticPressure(
                    tuple(loaded_plate_ids), direction_code, 2,
                    float(direction_sign * load_at_bottom), float(direction_sign * load_at_top)
                )
                if success:
                    applied_natively = True
            except Exception:
                pass

        if not applied_natively:
            print(f"Native command not confirmed for '{load_name}'. Applying precise discrete per-plate loads...")
            applied_count = 0
            for pid, y1, y2 in stem_plates_for_load:
                if y1 >= pressure_height:
                    continue
                y_top_eff = min(y2, pressure_height)
                y_avg = (y1 + y_top_eff) / 2.0
                depth = pressure_height - y_avg
                p = load_at_top + (load_at_bottom - load_at_top) * (depth / pressure_height)
                p_signed = direction_sign * p
                ok = False
                try:
                    ok = load_interface.AddElementPressure([int(pid)], direction_code, float(p_signed), 0, 0, 0, 0)
                except Exception:
                    try:
                        ok = load_interface.AddElementPressure(int(pid), direction_code, float(p_signed), 0, 0, 0, 0)
                    except Exception as e:
                        print(f"Failed to apply discrete pressure to plate {pid}: {e}")
                if ok:
                    applied_count += 1
            print(f"Discrete '{load_name}' pressure applied to {applied_count}/{len(stem_plates_for_load)} stem plates.")

    except Exception as e:
        QMessageBox.critical(None, "Pressure Load Error", f"A critical error occurred applying '{load_name}':\n{e}")


def apply_hydrostatic_load(staad, stem_plates_for_load, load_case, load_at_top, load_at_bottom, pressure_height,
                           direction_code=3, direction_sign=-1, load_name="Water Hydrostatic"):
    return apply_triangular_pressure_load(
        staad, stem_plates_for_load, load_case, load_at_top, load_at_bottom, pressure_height,
        direction_code=direction_code, direction_sign=direction_sign, load_name=load_name
    )


# [NEW] ------------------------------------------------------------
#  SOIL SURCHARGE on the base slab, defined separately for HEEL and TOE.
#  It is added INSIDE the dead load case (load_case = LC1, the self-weight
#  case), not as a separate load case.
#  Entered as a pressure in kN/m2, POSITIVE = DOWNWARD (gravity).
#  Applied over the full plate area in global Y, so the value written to
#  STAAD is -pressure. A negative entry acts upward.
#  If the dead load case does not exist yet (self-weight unticked) it is
#  created here. Returns the number of plates loaded; nothing is created
#  when both pressures are 0 or the plate lists are empty.
# ------------------------------------------------------------------
def apply_surcharge_load(staad, heel_plate_ids, toe_plate_ids, load_case,
                         heel_pressure, toe_pressure, load_name="Self Weight + Surcharge"):
    jobs = []
    for label, ids, p in (("heel", heel_plate_ids, heel_pressure),
                          ("toe", toe_plate_ids, toe_pressure)):
        plate_ids = [int(i) for i in ids]
        if plate_ids and abs(float(p)) > 1e-12:
            jobs.append((label, plate_ids, float(p)))
    if not jobs:
        print("Soil surcharge: nothing to apply (zero pressure or no plates in that zone).")
        return 0

    total = 0
    try:
        load_interface = staad.Load
        active_ok = False
        try:
            active_ok = bool(load_interface.SetLoadActive(int(load_case)))
        except Exception:
            active_ok = False
        if not active_ok:  # dead load case not created yet (self-weight off)
            try:
                load_interface.CreateNewPrimaryLoadEx2(str(load_name), 1, int(load_case))
            except Exception:
                pass
            load_interface.SetLoadActive(int(load_case))

        for label, plate_ids, p in jobs:
            value = -p  # downward in global Y
            batch_ok = False
            for arg in (plate_ids, tuple(plate_ids)):
                try:
                    batch_ok = bool(load_interface.AddElementPressure(
                        arg, SURCHARGE_DIRECTION_CODE, float(value), 0, 0, 0, 0))
                except Exception:
                    batch_ok = False
                if batch_ok:
                    break

            done = 0
            if batch_ok:
                done = len(plate_ids)
            else:
                for pid in plate_ids:
                    ok = False
                    for arg in ([pid], pid):
                        try:
                            ok = bool(load_interface.AddElementPressure(
                                arg, SURCHARGE_DIRECTION_CODE, float(value), 0, 0, 0, 0))
                        except Exception:
                            ok = False
                        if ok:
                            break
                    if ok:
                        done += 1
            print(f"Soil surcharge on {label}: {p:g} kN/m2 (downward) added to "
                  f"{done}/{len(plate_ids)} plates in dead load case {load_case}.")
            total += done
    except Exception as e:
        QMessageBox.critical(None, "Surcharge Load Error",
                             f"A critical error occurred applying the soil surcharge:\n{e}")
    return total


# ==================================================================
#  MAIN STAAD MODEL GENERATOR
# ==================================================================
def generate_retaining_wall(mesh_type, mesh_size, length_z, toe_x, heel_width, height_y, counterfort_spacing,
                            t_base, t_stem_top, t_stem_bot, t_cf, apply_hydro, hydro_params,
                            curve_radius=None, curve_turn_dir=1,
                            apply_soil_support=True, soil_subgrade=DEFAULT_SOIL_SUBGRADE_KN_M3,
                            apply_selfweight=True, selfweight_load_case=SELFWEIGHT_LOAD_CASE,
                            apply_lateral_soil_support=True,
                            lateral_soil_subgrade=DEFAULT_LATERAL_SUBGRADE_KN_M3,
                            lateral_direction=DIRECTION_X,
                            lateral_support_mode=DEFAULT_LATERAL_SUPPORT_MODE,
                            lateral_compression_only=True,
                            curved_lateral_method=DEFAULT_CURVED_LATERAL_METHOD,
                            strut_length=DEFAULT_STRUT_LENGTH_M,
                            soil_side=DEFAULT_SOIL_SIDE,
                            apply_earth=False, earth_params=None,
                            repeat_load_specs=None,
                            counterfort_side_pattern=COUNTERFORT_SIDE_HEEL,
                            counterfort_mode=COUNTERFORT_MODE_UNIFORM,
                            counterfort_positions_str="",
                            perform_analysis_print_option=0,
                            counterfort_height=None,  # None / 0 = full wall height
                            save_folder=None,  # None = do not save to a chosen folder
                            save_filename=DEFAULT_SAVE_FILENAME,
                            save_method=DEFAULT_SAVE_METHOD,
                            open_after_save=False,
                            surcharge_params=None):  # [NEW] {'heel': kN/m2, 'toe': kN/m2}  (+ve = downward), added to LC1

    # Fast path: build the whole model in memory, write the .std file at the end
    direct_mode = bool(save_folder) and save_method == SAVE_METHOD_DIRECT

    if direct_mode:
        staad = StdModelRecorder(title="Retaining Wall")  # STAAD.Pro is not used while building
    else:
        if not os_analytical:
            QMessageBox.critical(None, "Library Error", "openstaadpy is not installed.")
            return
        try:
            staad = os_analytical.connect()
        except Exception as e:
            QMessageBox.critical(None, "Connection Error", f"Error connecting to STAAD.Pro:\n{e}")
            return

    stem_x = heel_width
    width_x = toe_x + heel_width

    use_curve = bool(curve_radius) and curve_radius > 0
    if use_curve:
        eff_r_heel = curve_radius - curve_turn_dir * (0 - stem_x)
        eff_r_toe = curve_radius - curve_turn_dir * (width_x - stem_x)
        if min(eff_r_toe, eff_r_heel) <= 0:
            QMessageBox.critical(None, "Curve Error", "Curve radius is too tight for the wall width.")
            return

    # Target file in the chosen folder
    saved_path = None
    save_dir = None
    if save_folder:
        try:
            os.makedirs(save_folder, exist_ok=True)
        except Exception as e:
            QMessageBox.critical(None, "Save Folder Error", f"Cannot use folder:\n{save_folder}\n\n{e}")
            return
        target_path = make_unique_path(build_model_path(save_folder, save_filename))
        saved_path = target_path
        save_dir = os.path.dirname(target_path)
        if not direct_mode:  # live method: create the file inside STAAD.Pro first
            print(f"Creating new STAAD file: {target_path}")
            staad, new_ok, new_err = create_new_staad_file(staad, target_path)
            if not new_ok:
                QMessageBox.critical(None, "Save Folder Error",
                                     f"Could not create the STAAD file:\n{target_path}\n\n{new_err}\n\n"
                                     "No geometry was generated.")
                return

    def curve_transform(x, y, z):
        if not use_curve:
            return x, y, z
        theta = z / curve_radius
        c, s = math.cos(theta), math.sin(theta)
        gx = stem_x + curve_turn_dir * curve_radius * (1 - c) + (x - stem_x) * c
        gz = s * (curve_radius - curve_turn_dir * (x - stem_x))
        return gx, y, gz

    node_dict = {}
    node_coords_global = {}

    def get_or_create_node(x, y, z):
        coord = (round(x, 4), round(y, 4), round(z, 4))
        if coord not in node_dict:
            gx, gy, gz = curve_transform(x, y, z)
            node_id = staad.Geometry.AddNode(gx, gy, gz)
            node_dict[coord] = node_id
            node_coords_global[node_id] = (gx, gy, gz)
        return node_dict[coord]

    base_plates, cf_plates, stem_plates_for_load = [], [], []
    stem_row_plates = {}
    stem_node_trib = {}
    stem_node_coords = {}

    def add_plate(n1, n2, n3, n4, collection):
        pid = staad.Geometry.AddPlate(n1, n2, n3, n4)
        if pid and (isinstance(pid, int) or str(pid).isdigit()):
            plate_id = int(pid)
            collection.append(plate_id)
            return plate_id
        return None

    num_x_heel = max(1, int(round(heel_width / mesh_size)))
    num_x_toe = max(1, int(round(toe_x / mesh_size))) if toe_x > 0 else 0
    step_x_heel = heel_width / num_x_heel if num_x_heel > 0 else 0
    step_x_toe = toe_x / num_x_toe if num_x_toe > 0 else 0
    num_y_stem = max(1, int(round(height_y / mesh_size)))
    step_y = height_y / num_y_stem
    num_z = max(1, int(round(length_z / mesh_size)))
    step_z = length_z / num_z

    # Effective counterfort height (snapped to the stem mesh rows)
    cf_h = snap_counterfort_height(counterfort_height, height_y, step_y)
    print(f"Counterfort height: {cf_h:g} m (wall height {height_y:g} m).")

    counterfort_specs = resolve_counterfort_specs(
        counterfort_mode, counterfort_positions_str, counterfort_side_pattern,
        length_z, counterfort_spacing, step_z
    )

    print("Generating base and stem geometry...")
    for i in range(num_x_heel):
        for j in range(num_z):
            x1, x2 = i * step_x_heel, (i + 1) * step_x_heel
            z1, z2 = j * step_z, (j + 1) * step_z
            add_plate(get_or_create_node(x1, 0, z1), get_or_create_node(x2, 0, z1),
                      get_or_create_node(x2, 0, z2), get_or_create_node(x1, 0, z2), base_plates)

    n_heel_base = len(base_plates)  # [NEW] heel plates are created first, toe plates after

    if num_x_toe > 0:
        for i in range(num_x_toe):
            for j in range(num_z):
                x1, x2 = heel_width + i * step_x_toe, heel_width + (i + 1) * step_x_toe
                z1, z2 = j * step_z, (j + 1) * step_z
                add_plate(get_or_create_node(x1, 0, z1), get_or_create_node(x2, 0, z1),
                          get_or_create_node(x2, 0, z2), get_or_create_node(x1, 0, z2), base_plates)

    # [NEW] base slab split into heel and toe plate sets (used by the soil surcharge)
    heel_base_plates = base_plates[:n_heel_base]
    toe_base_plates = base_plates[n_heel_base:]

    for j in range(num_z):
        for k in range(num_y_stem):
            y1, y2 = k * step_y, (k + 1) * step_y
            z1, z2 = j * step_z, (j + 1) * step_z
            t_avg = round(t_stem_bot - (t_stem_bot - t_stem_top) * ((y1 + y2) / 2.0 / height_y), 4)
            if t_avg not in stem_row_plates:
                stem_row_plates[t_avg] = []

            n1 = get_or_create_node(stem_x, y1, z1)
            n2 = get_or_create_node(stem_x, y1, z2)
            n3 = get_or_create_node(stem_x, y2, z2)
            n4 = get_or_create_node(stem_x, y2, z1)
            plate_id = add_plate(n1, n2, n3, n4, stem_row_plates[t_avg])
            if plate_id:
                stem_plates_for_load.append((plate_id, y1, y2))
                p1 = node_coords_global[n1]
                p2 = node_coords_global[n2]
                p3 = node_coords_global[n3]
                p4 = node_coords_global[n4]
                a_quarter = quad_area_3d(p1, p2, p3, p4) / 4.0
                for nn, pp in ((n1, p1), (n2, p2), (n3, p3), (n4, p4)):
                    stem_node_trib[nn] = stem_node_trib.get(nn, 0.0) + a_quarter
                    stem_node_coords[nn] = pp

    t_cf = float(t_cf)

    def get_side_geometry(side):
        # boundary lines use cf_h (counterfort height) instead of height_y
        if side == COUNTERFORT_SIDE_TOE:
            def boundary_f(x, y):
                return x * cf_h + y * toe_x - width_x * cf_h
            x_edges = [heel_width + i * step_x_toe for i in range(num_x_toe + 1)] if num_x_toe > 0 else []
            num_x_cells = num_x_toe
        else:
            def boundary_f(x, y):
                return heel_width * y - cf_h * x
            x_edges = [i * step_x_heel for i in range(num_x_heel + 1)]
            num_x_cells = num_x_heel

        def clip_rect_to_triangle(rect_pts):
            out = []
            n = len(rect_pts)
            for idx in range(n):
                curr = rect_pts[idx]
                nxt = rect_pts[(idx + 1) % n]
                f_curr = boundary_f(*curr)
                f_nxt = boundary_f(*nxt)
                curr_in = f_curr <= 1e-9
                next_in = f_nxt <= 1e-9
                if curr_in:
                    out.append(curr)
                if curr_in != next_in:
                    if abs(f_curr - f_nxt) > 1e-12:
                        t = f_curr / (f_curr - f_nxt)
                        ix = curr[0] + t * (nxt[0] - curr[0])
                        iy = curr[1] + t * (nxt[1] - curr[1])
                        out.append((round(ix, 4), round(iy, 4)))
            return clean_polygon(out)

        return boundary_f, x_edges, num_x_cells, clip_rect_to_triangle

    heel_geom = get_side_geometry(COUNTERFORT_SIDE_HEEL)
    toe_geom = get_side_geometry(COUNTERFORT_SIDE_TOE)

    def make_plates_from_polygon(poly_xy, zc, collection):
        poly = clean_polygon(poly_xy)
        if len(poly) < 3:
            return
        if len(poly) == 3:
            if triangle_area(poly[0], poly[1], poly[2]) > 1e-4:
                n1 = get_or_create_node(poly[0][0], poly[0][1], zc)
                n2 = get_or_create_node(poly[1][0], poly[1][1], zc)
                n3 = get_or_create_node(poly[2][0], poly[2][1], zc)
                add_plate(n1, n2, n3, 0, collection)
        elif len(poly) == 4:
            n1 = get_or_create_node(poly[0][0], poly[0][1], zc)
            n2 = get_or_create_node(poly[1][0], poly[1][1], zc)
            n3 = get_or_create_node(poly[2][0], poly[2][1], zc)
            n4 = get_or_create_node(poly[3][0], poly[3][1], zc)
            add_plate(n1, n2, n3, n4, collection)
        else:
            n0 = get_or_create_node(poly[0][0], poly[0][1], zc)
            for v in range(1, len(poly) - 1):
                if triangle_area(poly[0], poly[v], poly[v + 1]) > 1e-4:
                    nv = get_or_create_node(poly[v][0], poly[v][1], zc)
                    nv1 = get_or_create_node(poly[v + 1][0], poly[v + 1][1], zc)
                    add_plate(n0, nv, nv1, 0, collection)

    print("Generating counterforts...")
    for zc, side in counterfort_specs:
        zc = round(zc, 4)
        boundary_f, x_edges, num_x_cells, clip_rect_to_triangle = heel_geom if side == COUNTERFORT_SIDE_HEEL else toe_geom
        if num_x_cells <= 0:
            continue
        for i in range(num_x_cells):
            x1, x2 = x_edges[i], x_edges[i + 1]
            for k in range(num_y_stem):
                y_bot, y_top = k * step_y, (k + 1) * step_y
                rect = [(x1, y_bot), (x2, y_bot), (x2, y_top), (x1, y_top)]
                f_vals = [boundary_f(x, y) for (x, y) in rect]
                if all(f <= 1e-9 for f in f_vals):
                    n1 = get_or_create_node(x1, y_bot, zc)
                    n2 = get_or_create_node(x2, y_bot, zc)
                    n3 = get_or_create_node(x2, y_top, zc)
                    n4 = get_or_create_node(x1, y_top, zc)
                    add_plate(n1, n2, n3, n4, cf_plates)
                elif all(f >= -1e-9 for f in f_vals):
                    continue
                else:
                    clipped = clip_rect_to_triangle(rect)
                    if len(clipped) >= 3:
                        make_plates_from_polygon(clipped, zc, cf_plates)

    print("Assigning plate thicknesses & materials...")
    prop = staad.Property
    if base_plates:
        base_thick_id = prop.CreatePlateThicknessProperty([t_base, t_base, t_base, t_base])
        prop.AssignPlateThickness(base_plates, base_thick_id)

    for t_avg, plate_ids in stem_row_plates.items():
        if not plate_ids:
            continue
        stem_thick_id = prop.CreatePlateThicknessProperty([t_avg, t_avg, t_avg, t_avg])
        prop.AssignPlateThickness(plate_ids, stem_thick_id)

    if cf_plates:
        cf_thick_id = prop.CreatePlateThicknessProperty([t_cf, t_cf, t_cf, t_cf])
        prop.AssignPlateThickness(cf_plates, cf_thick_id)

    all_plates = list(base_plates) + list(cf_plates)
    for row_ids in stem_row_plates.values():
        all_plates.extend(row_ids)

    material_ok = ensure_concrete_material(prop, CONCRETE_MATERIAL_NAME)
    if all_plates and material_ok:
        prop.AssignMaterialToPlate(CONCRETE_MATERIAL_NAME, all_plates)

    if apply_soil_support:
        print("Assigning base soil support...")
        add_soil_support(staad, base_plates, soil_subgrade)

    # --- Lateral soil ---
    used_curved_springs = False
    curved_spring_node_ids = []
    used_curved_members = False
    n_curved_members = 0

    if apply_lateral_soil_support:
        all_stem_plates = []
        for row_ids in stem_row_plates.values():
            all_stem_plates.extend(row_ids)

        use_members = (lateral_support_mode == LATERAL_MODE_AUTO_NORMAL) and (curved_lateral_method == CURVED_METHOD_MEMBERS)
        use_springs = (lateral_support_mode == LATERAL_MODE_AUTO_NORMAL) and (curved_lateral_method == CURVED_METHOD_SPRINGS) and use_curve

        if use_members:
            print(f"Assigning lateral soil support (compression-only members, soil_side={soil_side})...")
            n_curved_members, _ = add_lateral_soil_support_curved_members(
                staad, stem_node_trib, stem_node_coords,
                curve_radius, curve_turn_dir, stem_x, lateral_soil_subgrade,
                strut_length=strut_length, strut_E=CONCRETE_E_KN_M2,
                compression_only=lateral_compression_only,
                strut_material=CONCRETE_MATERIAL_NAME,
                soil_side=soil_side
            )
            used_curved_members = n_curved_members > 0
            if not used_curved_members:
                print("  -> No struts created; falling back to global Plate Mat.")
                add_lateral_soil_support(staad, all_stem_plates, lateral_soil_subgrade,
                                         lateral_direction, compression_only=lateral_compression_only)
        elif use_springs:
            print("Assigning lateral soil support (curved: inclined normal springs)...")
            n_applied, curved_spring_node_ids = add_lateral_soil_support_curved_normal(
                staad, stem_node_trib, stem_node_coords,
                curve_radius, curve_turn_dir, stem_x, lateral_soil_subgrade
            )
            used_curved_springs = n_applied > 0
            if n_applied == 0:
                print("  -> No inclined springs created; falling back to global Plate Mat.")
                add_lateral_soil_support(staad, all_stem_plates, lateral_soil_subgrade,
                                         lateral_direction, compression_only=lateral_compression_only)
        else:
            print("Assigning lateral soil support (global-X Plate Mat)...")
            add_lateral_soil_support(staad, all_stem_plates, lateral_soil_subgrade,
                                     lateral_direction, compression_only=lateral_compression_only)

    def _snap_pressure_height(requested):
        ph = requested if requested and requested > 0 else height_y
        ph = min(ph, height_y)
        snapped = round(round(ph / step_y) * step_y, 6)
        return max(step_y, min(snapped, height_y))

    # LC1 dead load (self-weight + soil surcharge), LC2 water, LC3 earth
    # (defined first so REPEATs can factor them)
    sur_heel_val = 0.0
    sur_toe_val = 0.0
    if surcharge_params:  # [NEW]
        sur_heel_val = float(surcharge_params.get('heel', 0.0))
        sur_toe_val = float(surcharge_params.get('toe', 0.0))
    surcharge_wanted = ((len(heel_base_plates) > 0 and abs(sur_heel_val) > 1e-12)
                        or (len(toe_base_plates) > 0 and abs(sur_toe_val) > 1e-12))
    dead_title = "Self Weight + Surcharge" if surcharge_wanted else "Self Weight"

    if apply_selfweight:
        print("Applying self-weight...")
        add_selfweight_load(staad, selfweight_load_case, title=dead_title)

    # [NEW] Soil surcharge on the base slab (heel / toe), added to the SAME dead load case
    surcharge_applied = False
    if surcharge_wanted:
        n_sur = apply_surcharge_load(
            staad, heel_base_plates, toe_base_plates, selfweight_load_case,
            sur_heel_val, sur_toe_val, load_name=dead_title
        )
        surcharge_applied = n_sur > 0

    if apply_hydro and hydro_params:
        snapped_height = _snap_pressure_height(hydro_params.get('height'))
        apply_triangular_pressure_load(
            staad=staad,
            stem_plates_for_load=stem_plates_for_load,
            load_case=hydro_params['load_case'],
            load_at_top=hydro_params['load_at_top'],
            load_at_bottom=hydro_params['load_at_bottom'],
            pressure_height=snapped_height,
            direction_code=hydro_params.get('direction_code', HYDRO_DIRECTION_OPTIONS[DEFAULT_HYDRO_DIRECTION]['code']),
            direction_sign=hydro_params.get('direction_sign', HYDRO_DIRECTION_OPTIONS[DEFAULT_HYDRO_DIRECTION]['sign']),
            load_name="Water Hydrostatic"
        )

    if apply_earth and earth_params:
        snapped_height_e = _snap_pressure_height(earth_params.get('height'))
        apply_triangular_pressure_load(
            staad=staad,
            stem_plates_for_load=stem_plates_for_load,
            load_case=earth_params['load_case'],
            load_at_top=earth_params['load_at_top'],
            load_at_bottom=earth_params['load_at_bottom'],
            pressure_height=snapped_height_e,
            direction_code=earth_params.get('direction_code', HYDRO_DIRECTION_OPTIONS[DEFAULT_HYDRO_DIRECTION]['code']),
            direction_sign=earth_params.get('direction_sign', HYDRO_DIRECTION_OPTIONS[DEFAULT_HYDRO_DIRECTION]['sign']),
            load_name="Lateral Earth Pressure"
        )

    # [NEW] The dead load case (LC1) exists only if self-weight and/or a surcharge was applied;
    # otherwise drop it from the REPEAT combinations
    dead_case_exists = bool(apply_selfweight) or surcharge_applied
    if repeat_load_specs and not dead_case_exists:
        kept_specs = []
        for spec in repeat_load_specs:
            kept_pairs = [p for p in spec['pairs'] if int(p[0]) != int(selfweight_load_case)]
            if kept_pairs:
                new_spec = dict(spec)
                new_spec['pairs'] = kept_pairs
                kept_specs.append(new_spec)
        repeat_load_specs = kept_specs

    # --- REPEAT LOAD design cases (created AFTER LC1/LC2/LC3 exist) ---
    if repeat_load_specs:
        print("Creating REPEAT LOAD design cases...")
        for spec in repeat_load_specs:
            add_repeat_load_case(
                staad,
                lc_num=spec['lc_num'],
                title=spec.get('title', f"Repeat LC{spec['lc_num']}"),
                load_factor_pairs=spec['pairs'],
                load_type=spec.get('load_type', LOAD_TYPE_DEAD)
            )

    perform_analysis_ok = False
    try:
        staad.Command.PerformAnalysis(int(perform_analysis_print_option))
        perform_analysis_ok = True
    except Exception as e:
        QMessageBox.critical(None, "PERFORM ANALYSIS Error", f"Could not add the PERFORM ANALYSIS command:\n{e}")

    # --- Curved inclined springs, direct mode: SPRING COMPRESSION goes straight into the file ---
    if direct_mode and used_curved_springs and lateral_compression_only and curved_spring_node_ids:
        staad.AddSpringCompression(curved_spring_node_ids)

    # Save the model into the chosen folder
    save_msg = ""
    if saved_path and direct_mode:  # write the in-memory model as one .std file
        try:
            staad.write(saved_path)
            print(f"Model written to: {saved_path}")
            save_msg = f"\n\nModel written to:\n{saved_path}\n({staad.summary()})"
        except Exception as e:
            save_msg = f"\n\nCould not write the .std file:\n{saved_path}\n{e}"
        else:
            if open_after_save:
                open_ok, open_err = open_std_in_staad(saved_path)
                if open_ok:
                    save_msg += "\nOpened in STAAD.Pro."
                else:
                    save_msg += f"\nCould not open it in STAAD.Pro ({open_err}). Open the file manually."
    elif saved_path:
        save_ok, save_err = save_staad_model(staad)
        if save_ok:
            print(f"Model saved to: {saved_path}")
            save_msg = f"\n\nModel saved to:\n{saved_path}"
        else:
            print(f"SaveModel failed: {save_err}")
            save_msg = (f"\n\nModel file created at:\n{saved_path}\n"
                        f"but SaveModel failed ({save_err}). Press Ctrl+S in STAAD.Pro to save it.")

    extra_msg = ""
    if repeat_load_specs:
        lines = []
        for spec in repeat_load_specs:
            desc = " + ".join(f"{f:g}xLC{lc}" for lc, f in spec['pairs'])
            lines.append(f"  LC{spec['lc_num']} '{spec.get('title','')}' = {desc}")
        extra_msg += "\n\nREPEAT LOAD design cases created:\n" + "\n".join(lines)

    if surcharge_applied:  # [NEW]
        extra_msg += (f"\n\nSoil surcharge added to the dead load case (LC{selfweight_load_case}, "
                      f"positive = downward): heel {sur_heel_val:g} kN/m2, toe {sur_toe_val:g} kN/m2.")

    if used_curved_members and lateral_compression_only:
        extra_msg += (f"\n\nCompression-only soil struts (soil on {soil_side}): {n_curved_members} "
                      "MEMBER COMPRESSION members created, each with material + section + a FIXED anchor.")
    elif direct_mode and used_curved_springs and lateral_compression_only and curved_spring_node_ids:
        extra_msg += (f"\n\nCOMPRESSION-ONLY (curved springs): SET NL, SPRING COMPRESSION for "
                      f"{len(curved_spring_node_ids)} sprung nodes and CHANGE were written into the .std file.")
    elif used_curved_springs and lateral_compression_only and curved_spring_node_ids:
        n_primary = 0
        if apply_selfweight:
            n_primary += 1
        if apply_hydro and hydro_params:
            n_primary += 1
        if apply_earth and earth_params:
            n_primary += 1
        if repeat_load_specs:
            n_primary += len(repeat_load_specs)
        n_primary = max(1, n_primary)
        _, cmd_block, file_path = apply_or_emit_spring_compression(
            staad, curved_spring_node_ids, n_load_cases=n_primary,
            tc_mode=SPRING_TC_MODE, direction_token=SPRING_TC_DIRECTION_TOKEN,
            output_dir=(save_dir or os.getcwd())
        )
        loc = f"\nSaved to:\n{file_path}" if file_path else ""
        extra_msg += ("\n\nCOMPRESSION-ONLY (curved springs):\nOpenSTAAD cannot set SPRING COMPRESSION "
                      f"directly, so the command block for {len(curved_spring_node_ids)} sprung nodes was "
                      "generated for pasting into the STAAD editor." + loc)
        box = QMessageBox(None)
        box.setIcon(QMessageBox.Information)
        box.setWindowTitle("SPRING COMPRESSION commands (paste into STAAD editor)")
        box.setText("OpenSTAAD cannot write SPRING COMPRESSION.\nPaste the block below into the STAAD editor.")
        box.setDetailedText(cmd_block)
        box.exec()

    extra_msg += save_msg

    if perform_analysis_ok:
        QMessageBox.information(None, "Success",
                               "Model generation process is complete. PERFORM ANALYSIS command included." + extra_msg)
    else:
        QMessageBox.warning(None, "Completed with Warning",
                            "Model generation is complete, but the PERFORM ANALYSIS command could NOT be added. "
                            "See the error message above." + extra_msg)


# ==================================================================
#  PYVISTA PREVIEW MESH BUILDER
# ==================================================================
def build_wall_mesh(height, toe_x, heel_width, length_z,
                    t_base, t_stem_top, t_stem_bot, t_cf,
                    cf_spacing, cf_side_pattern=COUNTERFORT_SIDE_HEEL,
                    curve_radius=None, curve_turn_dir=1,
                    apply_hydro=False, hydro_min=0.0, hydro_max=60.0, hydro_height=None, hydro_dir="Local Z (Normal to Wall)",
                    apply_earth=False, earth_min=0.0, earth_max=40.0, earth_height=None, earth_dir="Local Z (Normal to Wall)",
                    mesh_size=None, cf_mode=COUNTERFORT_MODE_UNIFORM, cf_positions_str="",
                    show_struts=False, strut_length=DEFAULT_STRUT_LENGTH_M,
                    lateral_support_mode=DEFAULT_LATERAL_SUPPORT_MODE,
                    curved_lateral_method=DEFAULT_CURVED_LATERAL_METHOD,
                    apply_lateral_soil=True, soil_side=DEFAULT_SOIL_SIDE,
                    cf_height=None):  # None / 0 = full wall height

    stem_x = heel_width
    width_x = toe_x + heel_width
    meshes = []

    num_z_segments = 50
    use_curve = bool(curve_radius) and curve_radius > 0

    def curve_transform_points(points_array):
        if not use_curve:
            return points_array
        transformed_points = []
        for p in points_array:
            x, y, z = p
            theta = z / curve_radius
            c, s = math.cos(theta), math.sin(theta)
            gx = stem_x + curve_turn_dir * curve_radius * (1 - c) + (x - stem_x) * c
            gz = s * (curve_radius - curve_turn_dir * (x - stem_x))
            transformed_points.append([gx, y, gz])
        return np.array(transformed_points)

    def local_point(x_local, y_local, z_local):
        if not use_curve:
            return np.array([x_local, y_local, z_local])
        theta = z_local / curve_radius
        c, s = math.cos(theta), math.sin(theta)
        gx = stem_x + curve_turn_dir * curve_radius * (1 - c) + (x_local - stem_x) * c
        gz = s * (curve_radius - curve_turn_dir * (x_local - stem_x))
        return np.array([gx, y_local, gz])

    def offset_point(x_local, y_local, z_local, offset, is_local_z, global_vec, side_sign=1.0):
        if is_local_z:
            return local_point(x_local - side_sign * offset, y_local, z_local)
        face_pt = local_point(x_local, y_local, z_local)
        return face_pt + side_sign * global_vec * offset

    def build_pressure_diagram(stem_x_local, h_p, length_z, p_min, p_max, dir_str,
                               factor=0.05, nz=60, side_sign=1.0):
        if h_p <= 0:
            return None
        is_local_z = "Local Z" in dir_str
        if "Global Y" in dir_str:
            global_vec = np.array([0.0, 1.0, 0.0])
        elif "Global Z" in dir_str:
            global_vec = np.array([0.0, 0.0, 1.0])
        else:
            global_vec = np.array([1.0, 0.0, 0.0])

        max_offset = p_max * factor
        min_offset = p_min * factor
        z_vals = np.linspace(0, length_z, nz)

        base_face, base_outer, top_face, top_outer = [], [], [], []
        for z in z_vals:
            base_face.append(local_point(stem_x_local, 0.0, z))
            base_outer.append(offset_point(stem_x_local, 0.0, z, max_offset, is_local_z, global_vec, side_sign))
            top_face.append(local_point(stem_x_local, h_p, z))
            top_outer.append(offset_point(stem_x_local, h_p, z, min_offset, is_local_z, global_vec, side_sign))

        def ribbon(line_a, line_b, p_a, p_b):
            pa, pb = np.array(line_a), np.array(line_b)
            X, Y, Z = np.vstack([pa[:, 0], pb[:, 0]]), np.vstack([pa[:, 1], pb[:, 1]]), np.vstack([pa[:, 2], pb[:, 2]])
            grid = pv.StructuredGrid(X, Y, Z)
            grid["pressure"] = np.vstack([np.full(len(pa), p_a), np.full(len(pb), p_b)]).ravel(order="F")
            return grid.extract_geometry()

        return pv.MultiBlock([
            ribbon(base_face, base_outer, p_max, p_max),
            ribbon(base_outer, top_outer, p_max, p_min),
            ribbon(top_outer, top_face, p_min, p_min)
        ]).combine(merge_points=True)

    side_sign = -1.0 if soil_side == SOIL_SIDE_LEFT else 1.0

    xxb, yyb, zzb = np.meshgrid(np.array([0, width_x]), np.array([-t_base, 0]), np.linspace(0, length_z, num_z_segments), indexing='ij')
    base_grid = pv.StructuredGrid(xxb, yyb, zzb)
    base_grid.points = curve_transform_points(base_grid.points)
    meshes.append(base_grid.extract_geometry())

    xxs, yys, zzs = np.meshgrid(np.array([stem_x - t_stem_bot, stem_x]), np.array([0, height]), np.linspace(0, length_z, num_z_segments), indexing='ij')
    stem_grid = pv.StructuredGrid(xxs, yys, zzs)
    stem_grid.points = curve_transform_points(stem_grid.points)
    meshes.append(stem_grid.extract_geometry())

    step_z_preview = None
    step_y_preview = None
    if mesh_size and mesh_size > 0:
        num_z_mesh = max(1, int(round(length_z / mesh_size)))
        step_z_preview = length_z / num_z_mesh
        num_y_mesh_prev = max(1, int(round(height / mesh_size)))
        step_y_preview = height / num_y_mesh_prev

    # Effective counterfort height (same snapping as the generator)
    cf_h = snap_counterfort_height(cf_height, height, step_y_preview)

    cf_specs = resolve_counterfort_specs(
        cf_mode, cf_positions_str, cf_side_pattern, length_z, cf_spacing, step_z_preview
    )

    for zc, side in cf_specs:
        if side == COUNTERFORT_SIDE_TOE:
            tri_pts = np.array([[stem_x, 0, zc], [width_x, 0, zc], [stem_x, cf_h, zc]])
        else:
            tri_pts = np.array([[stem_x, 0, zc], [0, 0, zc], [stem_x, cf_h, zc]])
        gusset = pv.PolyData(tri_pts, np.hstack([[3, 0, 1, 2]])).extrude((0, 0, t_cf), capping=True)
        gusset.points = curve_transform_points(gusset.points)
        meshes.append(gusset)

    hydro_diagram = None
    if apply_hydro:
        h_p = min(hydro_height if hydro_height and hydro_height > 0 else height, height)
        hydro_diagram = build_pressure_diagram(stem_x, h_p, length_z, hydro_min, hydro_max,
                                               hydro_dir, factor=0.05, nz=60, side_sign=side_sign)

    earth_diagram = None
    if apply_earth:
        h_pe = min(earth_height if earth_height and earth_height > 0 else height, height)
        earth_diagram = build_pressure_diagram(stem_x, h_pe, length_z, earth_min, earth_max,
                                               earth_dir, factor=0.09, nz=60, side_sign=side_sign)

    strut_lines = None
    struts_applicable = (
        show_struts and apply_lateral_soil
        and lateral_support_mode == LATERAL_MODE_AUTO_NORMAL
        and curved_lateral_method == CURVED_METHOD_MEMBERS
    )
    if struts_applicable and mesh_size and mesh_size > 0:
        num_z_mesh = max(1, int(round(length_z / mesh_size)))
        step_z_m = length_z / num_z_mesh
        num_y_mesh = max(1, int(round(height / mesh_size)))
        step_y_m = height / num_y_mesh

        cx = stem_x + curve_turn_dir * curve_radius if use_curve else 0.0
        cz = 0.0

        pts = []
        lines = []
        idx = 0
        for jz in range(num_z_mesh + 1):
            zc = jz * step_z_m
            for ky in range(num_y_mesh + 1):
                yc = ky * step_y_m
                node = local_point(stem_x, yc, zc)
                gx, gy, gz = float(node[0]), float(node[1]), float(node[2])
                nx, nz, dist, points_inward = strut_direction(gx, gz, use_curve, cx, cz, soil_side)
                L = float(strut_length)
                if points_inward and L >= dist:
                    L = 0.5 * dist
                anchor = np.array([gx + nx * L, gy, gz + nz * L])
                pts.append([gx, gy, gz])
                pts.append([anchor[0], anchor[1], anchor[2]])
                lines.append([2, idx, idx + 1])
                idx += 2

        if pts:
            strut_lines = pv.PolyData(np.array(pts, dtype=float), lines=np.array(lines).ravel())

    combined = pv.MultiBlock(meshes).combine(merge_points=True)

    center_x = width_x / 2.0
    combined = combined.reflect((1, 0, 0), point=(center_x, 0, 0))
    if hydro_diagram is not None:
        hydro_diagram = hydro_diagram.reflect((1, 0, 0), point=(center_x, 0, 0))
    if earth_diagram is not None:
        earth_diagram = earth_diagram.reflect((1, 0, 0), point=(center_x, 0, 0))
    if strut_lines is not None:
        strut_lines = strut_lines.reflect((1, 0, 0), point=(center_x, 0, 0))

    return combined, hydro_diagram, earth_diagram, strut_lines


# ==================================================================
#  MAIN APPLICATION WINDOW
# ==================================================================
class WallGeneratorApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("STAAD.Pro Retaining Wall Generator -- PySide6 + PyVista")
        self.resize(1280, 800)

        splitter = QSplitter(Qt.Horizontal)

        self.tabs = QTabWidget()
        self.inputs = {}
        self.combos = {}
        self.checks = {}

        self._build_geometry_tab()
        self._build_thickness_tab()
        self._build_soil_tab()
        self._build_load_tab()

        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.addWidget(self.tabs)
        self._build_save_group(left_layout)

        btn_box = QGroupBox("Actions")
        btn_layout = QVBoxLayout(btn_box)
        row1 = QHBoxLayout()
        self.preview_btn = QPushButton("Update Preview")
        self.preview_btn.clicked.connect(self.update_preview)
        self.reset_btn = QPushButton("Reset View")
        self.reset_btn.clicked.connect(lambda: self.plotter.reset_camera())
        row1.addWidget(self.preview_btn)
        row1.addWidget(self.reset_btn)
        btn_layout.addLayout(row1)

        self.generate_btn = QPushButton("Generate STAAD Model")
        self.generate_btn.setStyleSheet("font-weight: bold; padding: 8px;")
        self.generate_btn.clicked.connect(self.on_generate_click)
        btn_layout.addWidget(self.generate_btn)

        left_layout.addWidget(btn_box)
        left_panel.setMaximumWidth(420)

        viewer_panel = QWidget()
        viewer_layout = QVBoxLayout(viewer_panel)
        self.plotter = QtInteractor(viewer_panel)
        viewer_layout.addWidget(self.plotter.interactor)

        splitter.addWidget(left_panel)
        splitter.addWidget(viewer_panel)
        splitter.setStretchFactor(1, 1)

        self.setCentralWidget(splitter)
        self.update_preview()

    def _build_save_group(self, parent_layout):
        box = QGroupBox("Save Model")
        form = QFormLayout(box)

        self.checks["save_model"] = QCheckBox("Save model to chosen folder")
        self.checks["save_model"].setChecked(False)
        self.checks["save_model"].toggled.connect(self._refresh_save_controls)
        form.addRow(self.checks["save_model"])

        row = QHBoxLayout()
        self.inputs["save_folder"] = QLineEdit("")
        self.inputs["save_folder"].setPlaceholderText("Choose a folder...")
        self.save_browse_btn = QPushButton("Browse...")
        self.save_browse_btn.clicked.connect(self.browse_save_folder)
        row.addWidget(self.inputs["save_folder"])
        row.addWidget(self.save_browse_btn)
        form.addRow("Folder:", row)

        self.inputs["save_name"] = QLineEdit(DEFAULT_SAVE_FILENAME)
        form.addRow("File name:", self.inputs["save_name"])

        # save method + open-after option
        self.combos["save_method"] = QComboBox()
        self.combos["save_method"].addItems(SAVE_METHOD_OPTIONS)
        self.combos["save_method"].setCurrentText(DEFAULT_SAVE_METHOD)
        self.combos["save_method"].currentTextChanged.connect(self._refresh_save_controls)
        form.addRow("Method:", self.combos["save_method"])

        self.checks["open_after"] = QCheckBox("Open saved model in STAAD.Pro")
        self.checks["open_after"].setChecked(True)
        form.addRow(self.checks["open_after"])

        parent_layout.addWidget(box)
        self._refresh_save_controls()

    def _refresh_save_controls(self, *args):
        on = self.checks["save_model"].isChecked()
        direct = self.combos["save_method"].currentText() == SAVE_METHOD_DIRECT
        self.inputs["save_folder"].setEnabled(on)
        self.inputs["save_name"].setEnabled(on)
        self.save_browse_btn.setEnabled(on)
        self.combos["save_method"].setEnabled(on)
        self.checks["open_after"].setEnabled(on and direct)

    def browse_save_folder(self):
        start = self.inputs["save_folder"].text().strip() or os.getcwd()
        folder = QFileDialog.getExistingDirectory(self, "Select folder to save the STAAD model", start)
        if folder:
            self.inputs["save_folder"].setText(os.path.normpath(folder))

    def _build_geometry_tab(self):
        tab = QWidget()
        form = QFormLayout(tab)
        fields = [
            ("mesh_size",   "Mesh Size (m):",       "0.5"),
            ("length_z",    "Wall Length (Z):",     "12.0"),
            ("toe_x",       "Toe Length (X):",      "1.5"),
            ("heel_width",  "Heel Width (X):",      "4.5"),
            ("height",      "Wall Height (Y):",     "6.0"),
            ("cf_spacing",  "Counterfort Spacing:", "4.0"),
            ("cf_height",   "Counterfort Height (0=Full):", "0.0"),
            ("curve_radius", "Curve Radius (0=No):", "0.0")
        ]
        for key, label, default in fields:
            edit = QLineEdit(default)
            self.inputs[key] = edit
            form.addRow(QLabel(label), edit)
        self.inputs["cf_height"].setPlaceholderText("0 = full wall height")

        self.combos["curve_dir"] = QComboBox()
        self.combos["curve_dir"].addItems(["Right (+X)", "Left (-X)"])
        form.addRow("Curve Direction:", self.combos["curve_dir"])

        self.combos["cf_side"] = QComboBox()
        self.combos["cf_side"].addItems(COUNTERFORT_SIDE_PATTERN_OPTIONS)
        form.addRow("Counterfort Side Pattern:", self.combos["cf_side"])

        self.combos["cf_mode"] = QComboBox()
        self.combos["cf_mode"].addItems([COUNTERFORT_MODE_UNIFORM, COUNTERFORT_MODE_MANUAL])
        form.addRow("Counterfort Placement:", self.combos["cf_mode"])

        self.inputs["cf_positions"] = QLineEdit("")
        self.inputs["cf_positions"].setPlaceholderText(
            "e.g. 3.0:heel, 6.5:toe, 9.0  (side optional; falls back to pattern above)"
        )
        form.addRow("Manual Positions (Z, m):", self.inputs["cf_positions"])

        self.tabs.addTab(tab, "Geometry")

    def _build_thickness_tab(self):
        tab = QWidget()
        form = QFormLayout(tab)
        fields = [
            ("t_base",     "Base Thickness (m):",     "0.5"),
            ("t_stem_bot", "Stem Bot Thickness (m):", "0.5"),
            ("t_stem_top", "Stem Top Thickness (m):", "0.3"),
            ("t_cf",       "Counterfort Thick (m):",  "0.4"),
        ]
        for key, label, default in fields:
            edit = QLineEdit(default)
            self.inputs[key] = edit
            form.addRow(QLabel(label), edit)
        self.tabs.addTab(tab, "Thickness")

    def _build_soil_tab(self):
        tab = QWidget()
        form = QFormLayout(tab)

        self.checks["soil"] = QCheckBox("Apply Base Soil Support")
        self.checks["soil"].setChecked(True)
        form.addRow(self.checks["soil"])
        self.inputs["soil_sub"] = QLineEdit(str(DEFAULT_SOIL_SUBGRADE_KN_M3))
        form.addRow("Base Subgrade (kN/m3):", self.inputs["soil_sub"])

        self.checks["lat_soil"] = QCheckBox("Apply Lateral Soil")
        self.checks["lat_soil"].setChecked(True)
        form.addRow(self.checks["lat_soil"])
        self.inputs["lat_sub"] = QLineEdit(str(DEFAULT_LATERAL_SUBGRADE_KN_M3))
        form.addRow("Lat. Subgrade (kN/m3):", self.inputs["lat_sub"])

        self.combos["soil_side"] = QComboBox()
        self.combos["soil_side"].addItems(SOIL_SIDE_OPTIONS)
        self.combos["soil_side"].setCurrentText(DEFAULT_SOIL_SIDE)
        form.addRow("Soil Side:", self.combos["soil_side"])

        self.combos["lat_mode"] = QComboBox()
        self.combos["lat_mode"].addItems(LATERAL_SUPPORT_MODE_OPTIONS)
        self.combos["lat_mode"].setCurrentText(DEFAULT_LATERAL_SUPPORT_MODE)
        form.addRow("Lateral Support Mode:", self.combos["lat_mode"])

        self.combos["curved_method"] = QComboBox()
        self.combos["curved_method"].addItems(CURVED_LATERAL_METHOD_OPTIONS)
        self.combos["curved_method"].setCurrentText(DEFAULT_CURVED_LATERAL_METHOD)
        form.addRow("Curved Lateral Method:", self.combos["curved_method"])

        self.inputs["strut_length"] = QLineEdit(str(DEFAULT_STRUT_LENGTH_M))
        form.addRow("Soil Strut Length (m):", self.inputs["strut_length"])

        self.checks["show_struts"] = QCheckBox("Show Soil Struts in Preview")
        self.checks["show_struts"].setChecked(True)
        form.addRow(self.checks["show_struts"])

        self.checks["lat_comp_only"] = QCheckBox("Lateral soil COMPRESSION-only")
        self.checks["lat_comp_only"].setChecked(True)
        form.addRow(self.checks["lat_comp_only"])

        self.strut_calc_label = QLabel()
        self.strut_calc_label.setWordWrap(True)
        self.strut_calc_label.setStyleSheet(
            "QLabel { background:#f4f6f8; border:1px solid #c8ccd0; border-radius:4px;"
            " padding:6px; font-family:Consolas,'Courier New',monospace; font-size:11px; }"
        )
        form.addRow("Equivalent Section:", self.strut_calc_label)
        for key in ("lat_sub", "strut_length", "mesh_size"):
            if key in self.inputs:
                self.inputs[key].textChanged.connect(self.update_strut_calc)
        self.combos["curved_method"].currentTextChanged.connect(self.update_strut_calc)
        self.combos["soil_side"].currentTextChanged.connect(self.update_preview)
        self.update_strut_calc()

        self.tabs.addTab(tab, "Supports")

    def update_strut_calc(self):
        try:
            ks = float(self.inputs["lat_sub"].text())
            L = float(self.inputs["strut_length"].text())
            ms = float(self.inputs["mesh_size"].text())
        except (ValueError, KeyError):
            self.strut_calc_label.setText("Enter numeric Lat. Subgrade, Strut Length and Mesh Size.")
            return
        if self.combos["curved_method"].currentText() != CURVED_METHOD_MEMBERS:
            self.strut_calc_label.setText("Applies to the 'Compression-Only Members' method.")
            return
        if ks <= 0 or L <= 0 or ms <= 0:
            self.strut_calc_label.setText("Lat. Subgrade, Strut Length and Mesh Size must be > 0.")
            return
        E = CONCRETE_E_KN_M2
        A_trib = ms * ms
        K = ks * A_trib
        A = K * L / E
        d_mm = math.sqrt(4.0 * A / math.pi) * 1000.0
        self.strut_calc_label.setText(
            f"K = k_s x mesh^2 = {ks:,.0f} x {A_trib:.3g} = {K:,.1f} kN/m\n"
            f"A = K x L / E    = {A:.3e} m2   (L={L:g} m, E={E:,.0f})\n"
            f"d = sqrt(4A/pi)  = {d_mm:.2f} mm"
        )

    def _build_load_tab(self):
        tab = QWidget()
        form = QFormLayout(tab)

        self.checks["selfweight"] = QCheckBox("Apply Self-Weight (LC1)")
        self.checks["selfweight"].setChecked(True)
        form.addRow(self.checks["selfweight"])

        # ---------------- WATER hydrostatic (LC2) ----------------
        self.checks["hydro"] = QCheckBox("Apply Hydrostatic Load - Water (LC2)")
        self.checks["hydro"].setChecked(True)
        self.checks["hydro"].stateChanged.connect(self.update_preview)
        form.addRow(self.checks["hydro"])

        self.inputs["hydro_max"] = QLineEdit("60.0")
        self.inputs["hydro_max"].textChanged.connect(self.update_preview)
        form.addRow("Max Water Pressure (kN/m2):", self.inputs["hydro_max"])

        self.inputs["hydro_height"] = QLineEdit("0.0")
        self.inputs["hydro_height"].setPlaceholderText("0 = full wall height")
        self.inputs["hydro_height"].textChanged.connect(self.update_preview)
        form.addRow("Water Height (m, 0=Full):", self.inputs["hydro_height"])

        self.combos["hydro_dir"] = QComboBox()
        self.combos["hydro_dir"].addItems(list(HYDRO_DIRECTION_OPTIONS.keys()))
        self.combos["hydro_dir"].setCurrentText(DEFAULT_HYDRO_DIRECTION)
        self.combos["hydro_dir"].currentTextChanged.connect(self.update_preview)
        form.addRow("Water Direction:", self.combos["hydro_dir"])

        # ---------------- LATERAL EARTH pressure (LC3) ----------------
        self.checks["earth"] = QCheckBox("Apply Lateral Earth Pressure - Soil (LC3)")
        self.checks["earth"].setChecked(False)
        self.checks["earth"].stateChanged.connect(self.update_preview)
        form.addRow(self.checks["earth"])

        self.inputs["earth_ka"] = QLineEdit(str(DEFAULT_EARTH_KA))
        self.inputs["earth_ka"].textChanged.connect(self.update_preview)
        form.addRow("Earth Coeff. Ka:", self.inputs["earth_ka"])

        self.inputs["earth_gamma"] = QLineEdit(str(DEFAULT_EARTH_GAMMA_KN_M3))
        self.inputs["earth_gamma"].textChanged.connect(self.update_preview)
        form.addRow("Soil Unit Weight (kN/m3):", self.inputs["earth_gamma"])

        self.inputs["earth_max"] = QLineEdit("0.0")
        self.inputs["earth_max"].setPlaceholderText("0 = auto = Ka*gamma*height; negative = reverse direction")
        self.inputs["earth_max"].textChanged.connect(self.update_preview)
        form.addRow("Max Earth Pressure (kN/m2):", self.inputs["earth_max"])

        self.inputs["earth_height"] = QLineEdit("0.0")
        self.inputs["earth_height"].setPlaceholderText("0 = full wall height")
        self.inputs["earth_height"].textChanged.connect(self.update_preview)
        form.addRow("Earth Height (m, 0=Full):", self.inputs["earth_height"])

        self.combos["earth_dir"] = QComboBox()
        self.combos["earth_dir"].addItems(list(HYDRO_DIRECTION_OPTIONS.keys()))
        self.combos["earth_dir"].setCurrentText(DEFAULT_HYDRO_DIRECTION)
        self.combos["earth_dir"].currentTextChanged.connect(self.update_preview)
        form.addRow("Earth Direction:", self.combos["earth_dir"])

        # ---------------- [NEW] SOIL SURCHARGE on base slab, heel / toe (inside Dead Load LC1) ----------------
        self.checks["surcharge"] = QCheckBox("Apply Soil Surcharge on Base - Heel / Toe (in Dead Load, LC1)")
        self.checks["surcharge"].setChecked(False)
        form.addRow(self.checks["surcharge"])

        self.inputs["sur_heel"] = QLineEdit("0.0")
        self.inputs["sur_heel"].setPlaceholderText("kN/m2, +ve = downward, 0 = none")
        form.addRow("Heel Surcharge (kN/m2, +ve down):", self.inputs["sur_heel"])

        self.inputs["sur_toe"] = QLineEdit("0.0")
        self.inputs["sur_toe"].setPlaceholderText("kN/m2, +ve = downward, 0 = none")
        form.addRow("Toe Surcharge (kN/m2, +ve down):", self.inputs["sur_toe"])

        # ---------------- REPEAT LOAD design cases ----------------
        combo_box = QGroupBox("REPEAT LOAD Design Cases (primary, P-Delta safe)")
        combo_form = QFormLayout(combo_box)

        # Combo 1: SW + Earth  (LC4)
        self.checks["combo1"] = QCheckBox("LC4: Self-weight + Earth")
        self.checks["combo1"].setChecked(True)
        combo_form.addRow(self.checks["combo1"])
        self.inputs["c1_sw"] = QLineEdit("1.0")
        combo_form.addRow("  SW factor:", self.inputs["c1_sw"])
        self.inputs["c1_earth"] = QLineEdit("1.0")
        combo_form.addRow("  Earth factor:", self.inputs["c1_earth"])

        # Combo 2: SW + Earth + Water  (LC5)
        self.checks["combo2"] = QCheckBox("LC5: Self-weight + Earth + Water")
        self.checks["combo2"].setChecked(True)
        combo_form.addRow(self.checks["combo2"])
        self.inputs["c2_sw"] = QLineEdit("1.0")
        combo_form.addRow("  SW factor:", self.inputs["c2_sw"])
        self.inputs["c2_earth"] = QLineEdit("1.0")
        combo_form.addRow("  Earth factor:", self.inputs["c2_earth"])
        self.inputs["c2_water"] = QLineEdit("1.0")
        combo_form.addRow("  Water factor:", self.inputs["c2_water"])

        form.addRow(combo_box)

        self.combos["analysis_print"] = QComboBox()
        self.combos["analysis_print"].addItems(list(PERFORM_ANALYSIS_PRINT_OPTIONS.keys()))
        self.combos["analysis_print"].setCurrentText(DEFAULT_PERFORM_ANALYSIS_PRINT_OPTION)
        form.addRow("Analysis Print Option:", self.combos["analysis_print"])

        self.tabs.addTab(tab, "Loads")

    def _earth_base_and_height(self, wall_h):
        try:
            ka = float(self.inputs["earth_ka"].text())
        except ValueError:
            ka = DEFAULT_EARTH_KA
        try:
            gamma = float(self.inputs["earth_gamma"].text())
        except ValueError:
            gamma = DEFAULT_EARTH_GAMMA_KN_M3
        try:
            emax = float(self.inputs["earth_max"].text())
        except ValueError:
            emax = 0.0
        try:
            eh = float(self.inputs["earth_height"].text())
        except ValueError:
            eh = 0.0
        eh = eh if eh > 0 else wall_h
        # Any non-zero entry is used as typed, INCLUDING negative values, so a
        # negative Max Earth Pressure reverses the load direction. Only 0 = automatic.
        base = emax if abs(emax) > 1e-12 else (ka * gamma * eh)
        return base, eh

    def _float_or(self, key, default):
        try:
            return float(self.inputs[key].text())
        except (ValueError, KeyError):
            return default

    def update_preview(self):
        try:
            h = float(self.inputs["height"].text())
            t_x = float(self.inputs["toe_x"].text())
            hw = float(self.inputs["heel_width"].text())
            lz = float(self.inputs["length_z"].text())
            tb = float(self.inputs["t_base"].text())
            tst = float(self.inputs["t_stem_top"].text())
            tsb = float(self.inputs["t_stem_bot"].text())
            tcf = float(self.inputs["t_cf"].text())
            cfs = float(self.inputs["cf_spacing"].text())
            cr = float(self.inputs["curve_radius"].text())
            ms = float(self.inputs["mesh_size"].text())
            cfh = self._float_or("cf_height", 0.0)

            cdir_text = self.combos["curve_dir"].currentText()
            cdir = 1 if "Right" in cdir_text else -1

            cf_mode_sel = self.combos["cf_mode"].currentText()
            cf_positions_sel = self.inputs["cf_positions"].text()
            cf_side_pattern_sel = self.combos["cf_side"].currentText()

            apply_hyd = self.checks["hydro"].isChecked()
            hyd_max = float(self.inputs["hydro_max"].text())
            hyd_height_val = float(self.inputs["hydro_height"].text())
            hydro_dir_sel = self.combos["hydro_dir"].currentText()

            apply_earth_p = self.checks["earth"].isChecked()
            earth_base, earth_h = self._earth_base_and_height(h)
            earth_dir_sel = self.combos["earth_dir"].currentText()

            lat_mode_sel = self.combos["lat_mode"].currentText()
            curved_method_sel = self.combos["curved_method"].currentText()
            soil_side_sel = self.combos["soil_side"].currentText()
            apply_lat = self.checks["lat_soil"].isChecked()
            show_struts = self.checks["show_struts"].isChecked()
            try:
                strut_len = float(self.inputs["strut_length"].text())
            except ValueError:
                strut_len = DEFAULT_STRUT_LENGTH_M
            if strut_len <= 0:
                strut_len = DEFAULT_STRUT_LENGTH_M

            mesh, hydro_diag, earth_diag, strut_lines = build_wall_mesh(
                height=h, toe_x=t_x, heel_width=hw, length_z=lz,
                t_base=tb, t_stem_top=tst, t_stem_bot=tsb, t_cf=tcf,
                cf_spacing=cfs, cf_side_pattern=cf_side_pattern_sel,
                curve_radius=cr, curve_turn_dir=cdir,
                apply_hydro=apply_hyd, hydro_max=hyd_max,
                hydro_height=(hyd_height_val if hyd_height_val > 0 else None),
                hydro_dir=hydro_dir_sel,
                apply_earth=apply_earth_p, earth_max=earth_base,
                earth_height=(earth_h if earth_h > 0 else None),
                earth_dir=earth_dir_sel,
                mesh_size=ms, cf_mode=cf_mode_sel, cf_positions_str=cf_positions_sel,
                show_struts=show_struts, strut_length=strut_len,
                lateral_support_mode=lat_mode_sel,
                curved_lateral_method=curved_method_sel,
                apply_lateral_soil=apply_lat, soil_side=soil_side_sel,
                cf_height=cfh
            )

            self.plotter.clear()
            self.plotter.add_mesh(mesh, color="lightgray", show_edges=True)
            if hydro_diag is not None:
                self.plotter.add_mesh(hydro_diag, scalars="pressure", cmap="Blues",
                                      opacity=0.75, show_scalar_bar=True,
                                      scalar_bar_args={"title": "Water (kN/m2)"})
            if earth_diag is not None:
                self.plotter.add_mesh(earth_diag, scalars="pressure", cmap="Oranges",
                                      opacity=0.70, show_scalar_bar=True,
                                      scalar_bar_args={"title": "Earth (kN/m2)"})
            if strut_lines is not None and strut_lines.n_points > 0:
                self.plotter.add_mesh(strut_lines, color="red", line_width=3,
                                      render_lines_as_tubes=True)
                self.plotter.add_mesh(strut_lines.points, color="red",
                                      render_points_as_spheres=True, point_size=6)

            self.plotter.show_grid()
            self.plotter.reset_camera()

        except ValueError:
            QMessageBox.warning(self, "Input Error", "Please ensure all numeric inputs are valid.")

    def on_generate_click(self):
        try:
            ms = float(self.inputs["mesh_size"].text())
            lz = float(self.inputs["length_z"].text())
            t_x = float(self.inputs["toe_x"].text())
            hw = float(self.inputs["heel_width"].text())
            h = float(self.inputs["height"].text())
            cfs = float(self.inputs["cf_spacing"].text())
            cr = float(self.inputs["curve_radius"].text())
            cfh = self._float_or("cf_height", 0.0)

            # Save-to-folder option
            save_folder = None
            save_name = DEFAULT_SAVE_FILENAME
            save_method = DEFAULT_SAVE_METHOD
            open_after = False
            if self.checks["save_model"].isChecked():
                save_folder = self.inputs["save_folder"].text().strip()
                if not save_folder or not os.path.isdir(save_folder):
                    QMessageBox.warning(self, "Save Folder",
                                        "Please choose an existing folder to save the model "
                                        "(or untick 'Save model to chosen folder').")
                    return
                save_name = self.inputs["save_name"].text().strip() or DEFAULT_SAVE_FILENAME
                save_method = self.combos["save_method"].currentText()
                open_after = self.checks["open_after"].isChecked() and save_method == SAVE_METHOD_DIRECT

            cdir_text = self.combos["curve_dir"].currentText()
            cdir = 1 if "Right" in cdir_text else -1

            tb = float(self.inputs["t_base"].text())
            tst = float(self.inputs["t_stem_top"].text())
            tsb = float(self.inputs["t_stem_bot"].text())
            tcf = float(self.inputs["t_cf"].text())

            cf_mode_sel = self.combos["cf_mode"].currentText()
            cf_positions_sel = self.inputs["cf_positions"].text()
            cf_side_pattern_sel = self.combos["cf_side"].currentText()

            app_soil = self.checks["soil"].isChecked()
            soil_sub = float(self.inputs["soil_sub"].text())
            app_lat = self.checks["lat_soil"].isChecked()
            lat_sub = float(self.inputs["lat_sub"].text())
            lat_mode_sel = self.combos["lat_mode"].currentText()
            lat_comp_only = self.checks["lat_comp_only"].isChecked()
            curved_method_sel = self.combos["curved_method"].currentText()
            soil_side_sel = self.combos["soil_side"].currentText()
            try:
                strut_len = float(self.inputs["strut_length"].text())
            except ValueError:
                strut_len = DEFAULT_STRUT_LENGTH_M
            if strut_len <= 0:
                strut_len = DEFAULT_STRUT_LENGTH_M

            app_sw = self.checks["selfweight"].isChecked()
            app_hydro = self.checks["hydro"].isChecked()
            hyd_max = float(self.inputs["hydro_max"].text())
            hyd_height_val = float(self.inputs["hydro_height"].text())
            hydro_dir_sel = self.combos["hydro_dir"].currentText()
            hydro_dir_info = HYDRO_DIRECTION_OPTIONS.get(hydro_dir_sel, HYDRO_DIRECTION_OPTIONS[DEFAULT_HYDRO_DIRECTION])

            analysis_print_sel = self.combos["analysis_print"].currentText()
            analysis_print_option = PERFORM_ANALYSIS_PRINT_OPTIONS.get(
                analysis_print_sel, PERFORM_ANALYSIS_PRINT_OPTIONS[DEFAULT_PERFORM_ANALYSIS_PRINT_OPTION]
            )

            hydro_params = {
                'load_case': WATER_LOAD_CASE,
                'load_at_top': 0.0,
                'load_at_bottom': hyd_max,
                'height': (hyd_height_val if hyd_height_val > 0 else h),
                'direction_code': hydro_dir_info['code'],
                'direction_sign': hydro_dir_info['sign']
            }

            app_earth = self.checks["earth"].isChecked()
            earth_params = None
            if app_earth:
                earth_dir_sel = self.combos["earth_dir"].currentText()
                earth_dir_info = HYDRO_DIRECTION_OPTIONS.get(earth_dir_sel, HYDRO_DIRECTION_OPTIONS[DEFAULT_HYDRO_DIRECTION])
                earth_base, earth_h = self._earth_base_and_height(h)
                earth_params = {
                    'load_case': EARTH_LOAD_CASE,
                    'load_at_top': 0.0,
                    'load_at_bottom': earth_base,
                    'height': earth_h,
                    'direction_code': earth_dir_info['code'],
                    'direction_sign': earth_dir_info['sign']
                }

            # [NEW] Soil surcharge on the base slab, separate heel / toe pressures
            app_sur = self.checks["surcharge"].isChecked()
            surcharge_params = None
            if app_sur:
                surcharge_params = {
                    'heel': float(self.inputs["sur_heel"].text()),
                    'toe': float(self.inputs["sur_toe"].text())
                }

            # ---- Build REPEAT LOAD design cases from constituents that exist ----
            repeat_specs = []
            # Combo 1: SW (+ surcharge, both in LC1) + Earth
            if self.checks["combo1"].isChecked():
                pairs = []
                if app_sw or app_sur:  # [NEW] LC1 also carries the surcharge
                    pairs.append((SELFWEIGHT_LOAD_CASE, self._float_or("c1_sw", 1.0)))
                if app_earth:
                    pairs.append((EARTH_LOAD_CASE, self._float_or("c1_earth", 1.0)))
                if pairs:
                    repeat_specs.append({'lc_num': REPEAT_LC_SW_EARTH,
                                         'title': "SW + Earth", 'pairs': pairs})
            # Combo 2: SW (+ surcharge, both in LC1) + Earth + Water
            if self.checks["combo2"].isChecked():
                pairs = []
                if app_sw or app_sur:  # [NEW] LC1 also carries the surcharge
                    pairs.append((SELFWEIGHT_LOAD_CASE, self._float_or("c2_sw", 1.0)))
                if app_earth:
                    pairs.append((EARTH_LOAD_CASE, self._float_or("c2_earth", 1.0)))
                if app_hydro:
                    pairs.append((WATER_LOAD_CASE, self._float_or("c2_water", 1.0)))
                if pairs:
                    repeat_specs.append({'lc_num': REPEAT_LC_SW_EARTH_WATER,
                                         'title': "SW + Earth + Water", 'pairs': pairs})

            generate_retaining_wall(
                mesh_type="quad", mesh_size=ms, length_z=lz, toe_x=t_x, heel_width=hw,
                height_y=h, counterfort_spacing=cfs,
                t_base=tb, t_stem_top=tst, t_stem_bot=tsb, t_cf=tcf,
                apply_hydro=app_hydro, hydro_params=hydro_params,
                apply_earth=app_earth, earth_params=earth_params,
                repeat_load_specs=repeat_specs,
                curve_radius=cr, curve_turn_dir=cdir,
                apply_soil_support=app_soil, soil_subgrade=soil_sub,
                apply_selfweight=app_sw, selfweight_load_case=SELFWEIGHT_LOAD_CASE,
                apply_lateral_soil_support=app_lat, lateral_soil_subgrade=lat_sub,
                lateral_direction=DIRECTION_X,
                lateral_support_mode=lat_mode_sel,
                lateral_compression_only=lat_comp_only,
                curved_lateral_method=curved_method_sel,
                strut_length=strut_len,
                soil_side=soil_side_sel,
                counterfort_side_pattern=cf_side_pattern_sel,
                counterfort_mode=cf_mode_sel,
                counterfort_positions_str=cf_positions_sel,
                perform_analysis_print_option=analysis_print_option,
                counterfort_height=cfh,
                save_folder=save_folder,
                save_filename=save_name,
                save_method=save_method,
                open_after_save=open_after,
                surcharge_params=surcharge_params  # [NEW]
            )

        except ValueError:
            QMessageBox.warning(self, "Input Error", "Please ensure all numeric inputs are valid numbers.")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = WallGeneratorApp()
    window.show()
    sys.exit(app.exec())
