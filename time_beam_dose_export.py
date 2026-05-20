from connect import *
from datetime import datetime
import csv
import math
import os
import statistics
import time
import traceback
import numpy as np


# ============================================================
# Configuration
# ============================================================
N_REPEATS       = 5
CITRIX_USERNAME = "80030361"

# Candidate save destinations for M4. Paths that are not writable are
# skipped with a logged reason.
SAVE_PATHS = {
    "M4a_cwd"         : None,  # filled in at runtime with os.getcwd()
    "M4b_server_F"    : "F:/RayStationTiming_tmp",
    "M4c_citrix_C"    : "//Client/C$/Users/{user}/Documents/RayStationTiming_tmp".format(
                            user=CITRIX_USERNAME),
    "M4d_server_C"    : "C:/Temp/RayStationTiming_tmp",
}

# Test beam index (per spec)
TEST_BEAM_INDEX = 0


# ============================================================
# Helpers
# ============================================================
def stats(samples):
    """Return (mean, std) over a list of floats; std is sample stdev."""
    if not samples:
        return (float('nan'), float('nan'))
    if len(samples) == 1:
        return (samples[0], 0.0)
    return (statistics.mean(samples), statistics.stdev(samples))


def fmt_mean_std(mean_s, std_s):
    return "{:8.4f} +/- {:7.4f}".format(mean_s, std_s)


def try_make_dir(path):
    """Try to create a directory; return (ok, reason)."""
    try:
        os.makedirs(path, exist_ok=True)
        # Probe writability
        probe = os.path.join(path, ".write_probe_{}.tmp".format(os.getpid()))
        with open(probe, 'wb') as f:
            f.write(b"x")
        os.remove(probe)
        return (True, "ok")
    except Exception as exc:
        return (False, "unwritable: {}".format(exc))


def safe_remove(path):
    try:
        if path and os.path.exists(path):
            os.remove(path)
    except Exception:
        pass


def get_beam_dose(case, beam_set, scratch_log):
    """Return (bd, source_label).  Try FractionDose first; fall back to
    the first FractionEvaluation/DoseOnExamination/DoseEvaluation."""
    try:
        bd = beam_set.FractionDose.BeamDoses[TEST_BEAM_INDEX]
        # Touch a trivial attribute to confirm it is a live handle
        _ = bd.InDoseGrid
        return (bd, "FractionDose.BeamDoses[{}]".format(TEST_BEAM_INDEX))
    except Exception as exc_primary:
        scratch_log.append("  Primary BeamDose handle unavailable: {}".format(exc_primary))
        try:
            bd = (case.TreatmentDelivery
                      .FractionEvaluations[0]
                      .DoseOnExaminations[0]
                      .DoseEvaluations[0]
                      .BeamDoses[TEST_BEAM_INDEX])
            _ = bd.InDoseGrid
            label = ("case.TreatmentDelivery.FractionEvaluations[0]"
                     ".DoseOnExaminations[0].DoseEvaluations[0]"
                     ".BeamDoses[{}]".format(TEST_BEAM_INDEX))
            print("  NOTICE: Using fallback FractionEvaluation handle.")
            return (bd, label)
        except Exception as exc_fallback:
            raise RuntimeError(
                "Could not obtain a BeamDose handle.  Primary: {}  Fallback: {}"
                .format(exc_primary, exc_fallback)
            )


def fetch_grid_metadata(bd):
    """Pull nx/ny/nz, voxel size, corner from a BeamDose handle."""
    grid   = bd.InDoseGrid
    nx     = int(grid.NrVoxels.x)
    ny     = int(grid.NrVoxels.y)
    nz     = int(grid.NrVoxels.z)
    vx     = (float(grid.VoxelSize.x),
              float(grid.VoxelSize.y),
              float(grid.VoxelSize.z))
    corner = (float(grid.Corner.x),
              float(grid.Corner.y),
              float(grid.Corner.z))
    return nx, ny, nz, vx, corner


def time_block(fn, *args, **kwargs):
    """Run fn(*args, **kwargs), return (elapsed_seconds, result)."""
    t0 = time.perf_counter()
    result = fn(*args, **kwargs)
    t1 = time.perf_counter()
    return (t1 - t0, result)


def write_npy_durable(path, arr):
    """np.save + flush + fsync + close, to include durability in the timing."""
    with open(path, 'wb') as f:
        np.save(f, arr, allow_pickle=False)
        f.flush()
        os.fsync(f.fileno())


def write_npz_compressed(path, arr):
    np.savez_compressed(path, dose=arr)


# ============================================================
# Measurement runners
# ============================================================
def measure_m1_metadata(bd, n):
    times = []
    for _ in range(n):
        dt, _ = time_block(fetch_grid_metadata, bd)
        times.append(dt)
    return times


def measure_m2_handle_only(bd, n):
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        flat = bd.DoseValues.DoseData
        t1 = time.perf_counter()
        # Do NOT iterate.  Reference flat so it isn't optimised away.
        _ = flat
        times.append(t1 - t0)
    return times


def measure_m3a_list(bd, nz, ny, nx, n):
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        flat = bd.DoseValues.DoseData
        arr  = np.array(list(flat), dtype=np.float32).reshape(nz, ny, nx)
        t1 = time.perf_counter()
        times.append(t1 - t0)
        del arr
    return times


def measure_m3b_fromiter(bd, nz, ny, nx, n):
    times = []
    total = nz * ny * nx
    for _ in range(n):
        t0 = time.perf_counter()
        flat = bd.DoseValues.DoseData
        arr  = np.fromiter(flat, dtype=np.float32, count=total).reshape(nz, ny, nx)
        t1 = time.perf_counter()
        times.append(t1 - t0)
        del arr
    return times


def measure_m3c_tuple(bd, nz, ny, nx, n):
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        flat = bd.DoseValues.DoseData
        arr  = np.array(tuple(flat), dtype=np.float32).reshape(nz, ny, nx)
        t1 = time.perf_counter()
        times.append(t1 - t0)
        del arr
    return times


def measure_m3d_bulk(bd, nz, ny, nx, n):
    """Try a handful of bulk accessor names.  Return dict of
    {accessor_name: (times_list_or_None, note_str)}."""
    candidates = [
        ("DoseValues.DoseDataAsArray",   lambda: bd.DoseValues.DoseDataAsArray),
        ("DoseValues.GetDoseValues()",   lambda: bd.DoseValues.GetDoseValues()),
        ("GetDoseValues()",              lambda: bd.GetDoseValues()),
    ]
    results = {}
    for name, getter in candidates:
        # Probe once outside the timing loop
        try:
            probe = getter()
            _ = probe
        except AttributeError as exc:
            results[name] = (None, "not available in this build ({})".format(exc))
            continue
        except Exception as exc:
            results[name] = (None, "probe failed: {}".format(exc))
            continue

        times = []
        ok = True
        for _ in range(n):
            try:
                t0 = time.perf_counter()
                raw = getter()
                arr = np.asarray(raw, dtype=np.float32).reshape(nz, ny, nx)
                t1 = time.perf_counter()
                times.append(t1 - t0)
                del arr
            except Exception as exc:
                results[name] = (None, "runtime failure: {}".format(exc))
                ok = False
                break
        if ok:
            results[name] = (times, "ok")
    return results


def measure_m4_saves(dose_array, dest_paths, n, scratch_files):
    """For each writable destination, time .npy and .npz writes."""
    results = {}  # path_label -> {".npy": (times, size_bytes), ".npz": (...)}
    for label, path in dest_paths.items():
        ok, reason = try_make_dir(path)
        if not ok:
            results[label] = {"skipped": reason, "path": path}
            continue

        per_fmt = {"path": path}
        for ext, writer in (
            (".npy", write_npy_durable),
            (".npz", write_npz_compressed),
        ):
            times = []
            file_size = None
            target = os.path.join(path, "timing_dose{}".format(ext))
            try:
                for _ in range(n):
                    safe_remove(target)
                    t0 = time.perf_counter()
                    writer(target, dose_array)
                    t1 = time.perf_counter()
                    times.append(t1 - t0)
                    if file_size is None and os.path.exists(target):
                        file_size = os.path.getsize(target)
                    scratch_files.append(target)
                per_fmt[ext] = (times, file_size)
            except Exception as exc:
                per_fmt[ext] = (None, "write failed: {}".format(exc))
            finally:
                safe_remove(target)
        results[label] = per_fmt
    return results


def measure_m5_patient_save(patient, n):
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        patient.Save()
        t1 = time.perf_counter()
        times.append(t1 - t0)
    return times


def measure_m6_scriptable_export(case, beam_set, export_folder, n, scratch_files):
    """Time ScriptableDicomExport of the test beam's physical dose to
    export_folder.  Returns (times, note)."""
    beam_set_id = beam_set.BeamSetIdentifier()

    # First call: confirm the kwargs work in this build.
    pre_files = set(os.listdir(export_folder)) if os.path.isdir(export_folder) else set()
    try:
        t0 = time.perf_counter()
        case.ScriptableDicomExport(
            ExportFolderPath=export_folder,
            BeamSets=[beam_set_id],
            PhysicalBeamSetDoseForBeamSets=[beam_set_id],
            IgnorePreConditionWarnings=True,
        )
        t1 = time.perf_counter()
    except Exception as exc:
        return ([t1 - t0] if 't1' in dir() else [], "failed: {}".format(exc))

    times = [t1 - t0]
    # Track everything new for cleanup
    new_files = [os.path.join(export_folder, f)
                 for f in os.listdir(export_folder) if f not in pre_files]
    scratch_files.extend(new_files)

    for _ in range(n - 1):
        pre = set(os.listdir(export_folder))
        try:
            t0 = time.perf_counter()
            case.ScriptableDicomExport(
                ExportFolderPath=export_folder,
                BeamSets=[beam_set_id],
                PhysicalBeamSetDoseForBeamSets=[beam_set_id],
                IgnorePreConditionWarnings=True,
            )
            t1 = time.perf_counter()
            times.append(t1 - t0)
            new_now = [os.path.join(export_folder, f)
                       for f in os.listdir(export_folder) if f not in pre]
            scratch_files.extend(new_now)
        except Exception as exc:
            return (times, "failed on repeat: {}".format(exc))

    return (times, "ok")


# ============================================================
# Main
# ============================================================
scratch_files = []

try:
    SAVE_PATHS["M4a_cwd"] = os.getcwd()

    patient  = get_current("Patient")
    case     = get_current("Case")
    plan     = get_current("Plan")
    beam_set = get_current("BeamSet")

    try:
        rs_version = patient.RayStationVersion
    except Exception:
        rs_version = "unknown"

    print("=" * 70)
    print("RayStation per-beam dose readout/save timing")
    print("=" * 70)
    print("  Plan        : {}".format(plan.Name))
    print("  BeamSet id  : {}".format(beam_set.BeamSetIdentifier()))
    print("  Test beam   : {}".format(TEST_BEAM_INDEX))
    print("  N_REPEATS   : {}".format(N_REPEATS))
    print("  RS version  : {}".format(rs_version))
    print("  CWD         : {}".format(os.getcwd()))

    handle_notes = []
    bd, bd_source = get_beam_dose(case, beam_set, handle_notes)
    for note in handle_notes:
        print(note)
    print("  BeamDose src: {}".format(bd_source))

    # --- Metadata pulled ONCE, outside any timed readout block --------
    nx, ny, nz, voxel_size_cm, corner_cm = fetch_grid_metadata(bd)
    total_voxels = nx * ny * nz
    array_mb     = total_voxels * 4 / (1024.0 * 1024.0)
    print("  Grid        : {} x {} x {}  ({} voxels, {:.2f} MB float32)"
          .format(nx, ny, nz, total_voxels, array_mb))
    print("  Voxel (cm)  : x={:.4f} y={:.4f} z={:.4f}".format(*voxel_size_cm))
    print("  Corner (cm) : x={:.4f} y={:.4f} z={:.4f}".format(*corner_cm))
    print("")

    # ---- M1 --------------------------------------------------------
    print("[M1] Grid metadata fetch ...")
    m1_times = measure_m1_metadata(bd, N_REPEATS)

    # ---- M2 --------------------------------------------------------
    print("[M2] DoseData handle fetch (no iteration) ...")
    m2_times = measure_m2_handle_only(bd, N_REPEATS)

    # ---- M3 --------------------------------------------------------
    print("[M3a] Readout via np.array(list(flat)) ...")
    m3a_times = measure_m3a_list(bd, nz, ny, nx, N_REPEATS)

    print("[M3b] Readout via np.fromiter(flat, count=N) ...")
    m3b_times = measure_m3b_fromiter(bd, nz, ny, nx, N_REPEATS)

    print("[M3c] Readout via np.array(tuple(flat)) ...")
    m3c_times = measure_m3c_tuple(bd, nz, ny, nx, N_REPEATS)

    print("[M3d] Readout via bulk accessors ...")
    m3d_results = measure_m3d_bulk(bd, nz, ny, nx, N_REPEATS)

    # ---- Determine fastest readout for downstream M4/M6 ------------
    readout_candidates = {
        "M3a list"     : m3a_times,
        "M3b fromiter" : m3b_times,
        "M3c tuple"    : m3c_times,
    }
    for name, (times, note) in m3d_results.items():
        if times is not None:
            readout_candidates["M3d " + name] = times

    best_readout_name = None
    best_readout_mean = float('inf')
    for name, times in readout_candidates.items():
        m, _ = stats(times)
        if m == m or True:  # NaN guard not needed for non-empty lists
            if m < best_readout_mean:
                best_readout_mean = m
                best_readout_name = name

    print("  Fastest readout: {} ({:.4f} s)".format(best_readout_name, best_readout_mean))

    # Use M3b (fromiter) array as the canonical save subject if available;
    # otherwise fall back to a fresh M3a build.
    print("  Building reference dose array for save tests ...")
    flat = bd.DoseValues.DoseData
    try:
        dose_array = np.fromiter(flat, dtype=np.float32,
                                 count=total_voxels).reshape(nz, ny, nx)
    except Exception:
        flat = bd.DoseValues.DoseData
        dose_array = np.array(list(flat), dtype=np.float32).reshape(nz, ny, nx)

    # ---- M4 --------------------------------------------------------
    print("[M4] Save destinations (.npy durable + .npz compressed) ...")
    m4_results = measure_m4_saves(dose_array, SAVE_PATHS, N_REPEATS, scratch_files)

    # Identify best save destination (use .npy mean, fall back to .npz if .npy failed)
    best_save_label = None
    best_save_mean  = float('inf')
    best_save_path  = None
    best_save_ext   = None
    for label, info in m4_results.items():
        if "skipped" in info:
            continue
        for ext in (".npy", ".npz"):
            entry = info.get(ext)
            if not entry:
                continue
            times, _ = entry
            if times is None:
                continue
            m, _ = stats(times)
            if m < best_save_mean:
                best_save_mean  = m
                best_save_label = label
                best_save_path  = info["path"]
                best_save_ext   = ext

    if best_save_path is None:
        print("  No writable save destination succeeded.")
    else:
        print("  Fastest save: {} ({}) at {:.4f} s -> {}"
              .format(best_save_label, best_save_ext, best_save_mean, best_save_path))

    # ---- M5 --------------------------------------------------------
    print("[M5] patient.Save() overhead ...")
    m5_times = measure_m5_patient_save(patient, N_REPEATS)

    # ---- M6 --------------------------------------------------------
    if best_save_path is not None:
        print("[M6] ScriptableDicomExport to fastest save dir ...")
        m6_times, m6_note = measure_m6_scriptable_export(
            case, beam_set, best_save_path, N_REPEATS, scratch_files
        )
    else:
        m6_times = []
        m6_note  = "no writable destination"

    # ============================================================
    # Reporting
    # ============================================================
    print("")
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("  Plan          : {}".format(plan.Name))
    print("  BeamSet id    : {}".format(beam_set.BeamSetIdentifier()))
    print("  Test beam     : {}".format(TEST_BEAM_INDEX))
    print("  Grid          : {} x {} x {}  ({} voxels, {:.2f} MB float32)"
          .format(nx, ny, nz, total_voxels, array_mb))
    print("  Voxel (cm)    : x={:.4f} y={:.4f} z={:.4f}".format(*voxel_size_cm))
    print("  N_REPEATS     : {}".format(N_REPEATS))
    print("  RS version    : {}".format(rs_version))
    print("  BeamDose src  : {}".format(bd_source))
    print("")
    print("  {:<48s} | {:>22s} | Notes"
          .format("Measurement", "Mean +/- Std (s)"))
    print("  " + "-" * 96)

    def row(label, times, note=""):
        if times is None or len(times) == 0:
            print("  {:<48s} | {:>22s} | {}".format(label, "n/a", note))
        else:
            m, s = stats(times)
            print("  {:<48s} | {:>22s} | {}"
                  .format(label, fmt_mean_std(m, s), note))

    row("M1  Grid metadata fetch",                 m1_times)
    row("M2  DoseData handle fetch (no iteration)", m2_times)
    row("M3a Readout: np.array(list(flat))",       m3a_times)
    row("M3b Readout: np.fromiter(flat, count)",   m3b_times)
    row("M3c Readout: np.array(tuple(flat))",      m3c_times)
    for name, (times, note) in m3d_results.items():
        row("M3d Readout: " + name, times, note)

    for label, info in m4_results.items():
        if "skipped" in info:
            print("  {:<48s} | {:>22s} | {}  ({})"
                  .format("M4  " + label + "  (skipped)", "n/a",
                          info["skipped"], info["path"]))
            continue
        for ext in (".npy", ".npz"):
            entry = info.get(ext)
            if not entry:
                continue
            times, size = entry
            if times is None:
                row("M4  {} {}".format(label, ext), None, str(size))
            else:
                size_str = ("{:.2f} MB on disk".format(size / (1024.0 * 1024.0))
                            if size else "size unknown")
                row("M4  {} {}".format(label, ext), times,
                    "{}  ({})".format(size_str, info["path"]))

    row("M5  patient.Save()", m5_times)
    row("M6  ScriptableDicomExport (1 beam)", m6_times, m6_note)

    # ---- Footer analysis ------------------------------------------
    print("")
    print("=" * 70)
    print("FOOTER")
    print("=" * 70)

    m3a_mean, _ = stats(m3a_times)
    if best_readout_name is not None and m3a_mean > 0:
        speedup_readout = m3a_mean / best_readout_mean if best_readout_mean > 0 else float('inf')
        print("  Best readout      : {}  ({:.4f} s vs M3a {:.4f} s, {:.2f}x)"
              .format(best_readout_name, best_readout_mean, m3a_mean, speedup_readout))
    else:
        print("  Best readout      : unable to determine")

    citrix_info = m4_results.get("M4c_citrix_C", {})
    citrix_npz  = citrix_info.get(".npz") if isinstance(citrix_info, dict) else None
    citrix_mean = float('nan')
    if citrix_npz and citrix_npz[0]:
        citrix_mean, _ = stats(citrix_npz[0])

    if best_save_path is not None:
        if citrix_mean == citrix_mean and citrix_mean > 0 and best_save_mean > 0:
            sp_save = citrix_mean / best_save_mean
            print("  Best save dest    : {} ({})  {:.4f} s vs Citrix .npz {:.4f} s, {:.2f}x"
                  .format(best_save_label, best_save_ext,
                          best_save_mean, citrix_mean, sp_save))
        else:
            print("  Best save dest    : {} ({})  {:.4f} s  (Citrix baseline unavailable)"
                  .format(best_save_label, best_save_ext, best_save_mean))
    else:
        print("  Best save dest    : none succeeded")

    # Projected 6000-dose runs.  Current = M3a + Citrix .npz.
    current_per_dose = float('nan')
    if m3a_times and citrix_mean == citrix_mean:
        current_per_dose = m3a_mean + citrix_mean
    optimized_per_dose = float('nan')
    if best_readout_mean < float('inf') and best_save_mean < float('inf'):
        optimized_per_dose = best_readout_mean + best_save_mean

    def hours(per_dose_s, n=6000):
        if per_dose_s != per_dose_s:
            return float('nan')
        return per_dose_s * n / 3600.0

    cur_h = hours(current_per_dose)
    opt_h = hours(optimized_per_dose)
    if cur_h == cur_h and opt_h == opt_h and opt_h > 0:
        overall = cur_h / opt_h
        print("  6000-dose run     : current {:.2f} h  ->  optimized {:.2f} h  ({:.2f}x)"
              .format(cur_h, opt_h, overall))
    else:
        print("  6000-dose run     : current={}  optimized={}  (insufficient data)"
              .format(cur_h, opt_h))

    # ============================================================
    # CSV: one row per (measurement, repeat)
    # ============================================================
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path = os.path.join(os.getcwd(), "timing_results_{}.csv".format(ts))

    def csv_rows_for(label, times, notes=""):
        rows = []
        if times is None:
            rows.append([label, "n/a", "n/a", notes])
        else:
            for i, t in enumerate(times):
                rows.append([label, i, "{:.6f}".format(t), notes])
        return rows

    all_rows = []
    all_rows += csv_rows_for("M1_grid_metadata",   m1_times)
    all_rows += csv_rows_for("M2_handle_only",     m2_times)
    all_rows += csv_rows_for("M3a_list",           m3a_times)
    all_rows += csv_rows_for("M3b_fromiter",       m3b_times)
    all_rows += csv_rows_for("M3c_tuple",          m3c_times)
    for name, (times, note) in m3d_results.items():
        all_rows += csv_rows_for("M3d_" + name.replace(" ", "_"), times, note)
    for label, info in m4_results.items():
        if "skipped" in info:
            all_rows.append(["M4_" + label, "n/a", "n/a", info["skipped"]])
            continue
        for ext in (".npy", ".npz"):
            entry = info.get(ext)
            if not entry:
                continue
            times, size = entry
            note = "size_bytes={}".format(size) if isinstance(size, int) else str(size)
            all_rows += csv_rows_for("M4_{}{}".format(label, ext), times, note)
    all_rows += csv_rows_for("M5_patient_save",         m5_times)
    all_rows += csv_rows_for("M6_scriptable_export",    m6_times, m6_note)

    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["measurement", "repeat_index", "seconds", "notes"])
        # Header rows of context
        w.writerow(["#plan",         "",  "", plan.Name])
        w.writerow(["#beam_set_id",  "",  "", beam_set.BeamSetIdentifier()])
        w.writerow(["#test_beam",    "",  "", TEST_BEAM_INDEX])
        w.writerow(["#grid_nx_ny_nz","",  "", "{}x{}x{}".format(nx, ny, nz)])
        w.writerow(["#total_voxels", "",  "", total_voxels])
        w.writerow(["#voxel_cm",     "",  "", "{},{},{}".format(*voxel_size_cm)])
        w.writerow(["#n_repeats",    "",  "", N_REPEATS])
        w.writerow(["#rs_version",   "",  "", rs_version])
        w.writerow(["#bd_source",    "",  "", bd_source])
        for r in all_rows:
            w.writerow(r)

    print("")
    print("  CSV written: {}".format(csv_path))

except Exception as e:
    print("\nERROR: {}".format(str(e)))
    print(traceback.format_exc())

finally:
    # Clean up every scratch file we created.
    cleaned = 0
    for p in scratch_files:
        if p and os.path.exists(p):
            try:
                os.remove(p)
                cleaned += 1
            except Exception:
                pass
    if scratch_files:
        print("  Cleanup: removed {}/{} scratch file(s).".format(cleaned, len(scratch_files)))
