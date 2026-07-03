import csv
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import linear_sum_assignment


REPO = Path(__file__).resolve().parents[2]
CURRENT_MODEL = REPO / "results/mofa/mofa_res_10factor.hdf5"
SCALED_MODEL = REPO / "results/mofa/mofa_res_10factor_scale_views_true.hdf5"
OUTDIR = REPO / "results/mofa/scale_views_sensitivity"


def _decode(values):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def read_model(path):
    with h5py.File(path, "r") as handle:
        views = _decode(handle["views/views"][()])
        samples = _decode(handle["samples/single_group"][()])
        z = np.asarray(handle["expectations/Z/single_group"][()])
        r2_factor = np.asarray(handle["variance_explained/r2_per_factor/single_group"][()])
        r2_total = np.asarray(handle["variance_explained/r2_total/single_group"][()])
        elbo = np.asarray(handle["training_stats/elbo"][()])
        n_factors = np.asarray(handle["training_stats/number_factors"][()])

    return {
        "path": path,
        "views": views,
        "samples": samples,
        "z": z,
        "r2_factor": r2_factor,
        "r2_total": r2_total,
        "elbo": elbo,
        "n_factors": n_factors,
    }


def write_csv(path, rows, fieldnames):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def factor_correlations(current, scaled):
    if current["samples"] != scaled["samples"]:
        raise ValueError("Sample order differs between models.")

    corr = np.corrcoef(current["z"], scaled["z"])[: current["z"].shape[0], current["z"].shape[0] :]
    row_ind, col_ind = linear_sum_assignment(-np.abs(corr))
    rows = []
    for row, col in zip(row_ind, col_ind):
        rows.append(
            {
                "current_factor": f"Factor{row + 1}",
                "scale_views_true_factor": f"Factor{col + 1}",
                "pearson_r": corr[row, col],
                "abs_pearson_r": abs(corr[row, col]),
                "sign": 1 if corr[row, col] >= 0 else -1,
            }
        )
    rows.sort(key=lambda x: int(x["current_factor"].replace("Factor", "")))
    return rows, corr


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    current = read_model(CURRENT_MODEL)
    scaled = read_model(SCALED_MODEL)

    rows = []
    for label, model in [("current_scale_views_false", current), ("scale_views_true", scaled)]:
        for view, value in zip(model["views"], model["r2_total"]):
            rows.append({"model": label, "view": view, "r2_total_percent": value})
    write_csv(OUTDIR / "r2_total_by_view.csv", rows, ["model", "view", "r2_total_percent"])

    rows = []
    for label, model in [("current_scale_views_false", current), ("scale_views_true", scaled)]:
        for view_index, view in enumerate(model["views"]):
            for factor_index in range(model["r2_factor"].shape[1]):
                rows.append(
                    {
                        "model": label,
                        "view": view,
                        "factor": f"Factor{factor_index + 1}",
                        "r2_percent": model["r2_factor"][view_index, factor_index],
                    }
                )
    write_csv(OUTDIR / "r2_per_factor_by_view.csv", rows, ["model", "view", "factor", "r2_percent"])

    rows = []
    for label, model in [("current_scale_views_false", current), ("scale_views_true", scaled)]:
        for factor_index in range(model["r2_factor"].shape[1]):
            values = model["r2_factor"][:, factor_index]
            total = values.sum()
            dominant_index = int(np.argmax(values))
            rows.append(
                {
                    "model": label,
                    "factor": f"Factor{factor_index + 1}",
                    "dominant_view": model["views"][dominant_index],
                    "dominant_view_r2_percent": values[dominant_index],
                    "sum_r2_across_views_percent": total,
                    "dominant_view_fraction_of_factor_r2": values[dominant_index] / total if total else np.nan,
                }
            )
    write_csv(
        OUTDIR / "factor_view_dominance.csv",
        rows,
        [
            "model",
            "factor",
            "dominant_view",
            "dominant_view_r2_percent",
            "sum_r2_across_views_percent",
            "dominant_view_fraction_of_factor_r2",
        ],
    )

    match_rows, corr = factor_correlations(current, scaled)
    write_csv(
        OUTDIR / "factor_score_correlations.csv",
        match_rows,
        ["current_factor", "scale_views_true_factor", "pearson_r", "abs_pearson_r", "sign"],
    )
    np.savetxt(OUTDIR / "factor_score_correlation_matrix.csv", corr, delimiter=",")

    summary_rows = []
    for label, model in [("current_scale_views_false", current), ("scale_views_true", scaled)]:
        summary_rows.append(
            {
                "model": label,
                "path": str(model["path"].relative_to(REPO)),
                "active_factors": model["z"].shape[0],
                "iterations": len(model["elbo"]),
                "final_elbo": model["elbo"][-1],
                "final_recorded_factors": model["n_factors"][-1],
            }
        )
    write_csv(OUTDIR / "training_summary.csv", summary_rows, ["model", "path", "active_factors", "iterations", "final_elbo", "final_recorded_factors"])


if __name__ == "__main__":
    main()
