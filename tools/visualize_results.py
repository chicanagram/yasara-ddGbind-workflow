#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities to compare YASARA ddG outputs against experimental fold-change data.

The script is intentionally configured from the ``__main__`` block so it can be
run and tweaked easily from an IDE without adding CLI plumbing first.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils import drop_unnamed_columns, get_struct_variant_name, parse_csv_bool, parse_csv_list


DEFAULT_CALC_VALUE_COLUMNS = ["ebindDDG"]
DEFAULT_EXPERIMENTAL_VALUE_COLUMN = "fold_diff"
DEFAULT_EXPERIMENTAL_LOG_VALUE_COLUMN = "fold_diff_log10"
DEFAULT_CALC_ALIGN_COLUMNS = ["struct", "target_chain", "ligname", "mutations", "ff", "mvdist"]
DEFAULT_EXPERIMENTAL_ALIGN_COLUMNS = ["struct", "target_chain", "ligname", "mutations"]
DEFAULT_SEGMENT_COLUMNS = ["result_source", "ff", "mvdist"]
EXPERIMENTAL_METADATA_COLUMN_MAP = {
    "struct_name": "struct",
    "ligand_name": "ligname",
    "mutations": "mutations",
}
EXPERIMENTAL_SUPPORT_COLUMNS = [
    "pdb_id",
    "target_chain",
    "chain_id",
    "keep_multiple_chains_in_struct",
]
SUMMARY_SORT_COLUMNS = ["abs_spearman", "abs_pearson"]
CORE_OUTPUT_FILENAMES = {
    "calc_wide": "calculation_results_merged_wide.csv",
    "merged_chain_averaged": "calculation_vs_experiment_chain_averaged.csv",
    "summary_overall": "correlation_summary_overall.csv",
}
INTERMEDIATE_OUTPUT_FILENAMES = {
    "calc_long": "calculation_results_long.csv",
    "experimental": "experimental_values_normalized.csv",
    "merged_raw": "calculation_vs_experiment_raw.csv",
    "summary_raw": "correlation_summary_raw.csv",
    "summary_chain_averaged": "correlation_summary_chain_averaged.csv",
}
PLOT_MODES = {"grid", "individual"}


def normalize_string_columns(df: pd.DataFrame, columns) -> pd.DataFrame:
    df = df.copy()
    for column in columns:
        if column in df.columns:
            df[column] = df[column].fillna("").astype(str).str.strip()
    return df


def prepare_key_columns(df: pd.DataFrame, key_columns) -> pd.DataFrame:
    df = df.copy()
    for column in key_columns:
        if column not in df.columns:
            df[column] = pd.NA
    return normalize_string_columns(df, key_columns)


def infer_target_chain_from_struct(struct_value: object) -> str | None:
    if struct_value is None or pd.isna(struct_value):
        return None
    struct_text = str(struct_value).strip()
    if struct_text == "":
        return None

    match = re.search(r"_[^_]+-([A-Za-z0-9]+)_[^_]+$", struct_text)
    if match:
        return match.group(1)
    return None


def get_chain_agnostic_struct(struct_value: object, target_chain: object) -> str | None:
    if struct_value is None or pd.isna(struct_value):
        return None
    struct_text = str(struct_value).strip()
    if struct_text == "":
        return None

    chain_text = "" if target_chain is None or pd.isna(target_chain) else str(target_chain).strip()
    if chain_text == "":
        return struct_text

    pattern = rf"-{re.escape(chain_text)}_"
    return re.sub(pattern, "_", struct_text, count=1)


def check_unique_keys(df: pd.DataFrame, key_columns: list[str], label: str) -> None:
    duplicates = df[df.duplicated(subset=key_columns, keep=False)]
    if duplicates.empty:
        return

    preview = duplicates[key_columns].head(10).to_dict(orient="records")
    raise ValueError(
        f"{label} contains duplicate rows for alignment columns {key_columns}. "
        f"Adjust merge_id_columns or deduplicate upstream. Example duplicate keys: {preview}"
    )


def require_columns(df: pd.DataFrame, required_columns, label: str) -> None:
    missing_columns = [column for column in required_columns if column not in df.columns]
    if missing_columns:
        raise ValueError(f"{label} is missing required columns: {missing_columns}")


def load_calculation_file(
    file_path: Path,
    merge_id_columns: list[str],
    value_columns: list[str],
) -> pd.DataFrame:
    df = pd.read_csv(file_path)
    df = drop_unnamed_columns(df)
    df = df.drop_duplicates().reset_index(drop=True)

    if "target_chain" in merge_id_columns and "target_chain" not in df.columns and "struct" in df.columns:
        # Legacy files may encode the chain only in ``struct``.
        df["target_chain"] = df["struct"].map(infer_target_chain_from_struct)

    df = prepare_key_columns(df, merge_id_columns)
    if "struct" in df.columns:
        df["struct_base"] = [
            get_chain_agnostic_struct(struct_value, target_chain)
            for struct_value, target_chain in zip(df["struct"], df.get("target_chain", pd.Series(index=df.index)))
        ]
    require_columns(df, value_columns, str(file_path))

    keep_columns = list(dict.fromkeys(merge_id_columns + value_columns))
    df = df[keep_columns].copy()
    df["result_source"] = file_path.stem
    check_unique_keys(df, merge_id_columns, str(file_path))
    return df


def merge_calculation_results_wide(
    result_frames: list[pd.DataFrame],
    merge_id_columns: list[str],
    value_columns: list[str],
) -> pd.DataFrame:
    if not result_frames:
        return pd.DataFrame(columns=merge_id_columns)

    merged_df: pd.DataFrame | None = None
    for df in result_frames:
        source_label = df["result_source"].iloc[0]
        renamed = df.drop(columns=["result_source"]).rename(
            columns={column: f"{column}__{source_label}" for column in value_columns}
        )
        if merged_df is None:
            merged_df = renamed
        else:
            merged_df = merged_df.merge(renamed, on=merge_id_columns, how="outer", validate="one_to_one")

    assert merged_df is not None
    return merged_df.sort_values(merge_id_columns).reset_index(drop=True)


def concatenate_calculation_results_long(result_frames: list[pd.DataFrame]) -> pd.DataFrame:
    if not result_frames:
        return pd.DataFrame()
    return pd.concat(result_frames, axis=0, ignore_index=True, sort=False)


def load_experimental_data(
    experimental_file: Path,
    experimental_value_column: str,
    experimental_log_value_column: str,
    metadata_column_map: dict[str, str] | None = None,
) -> pd.DataFrame:
    metadata_column_map = metadata_column_map or EXPERIMENTAL_METADATA_COLUMN_MAP

    if not experimental_file.exists():
        raise FileNotFoundError(f"Experimental file not found: {experimental_file}")

    df = pd.read_csv(experimental_file)
    df = drop_unnamed_columns(df).rename(columns=metadata_column_map)
    required_columns = list(metadata_column_map.values()) + EXPERIMENTAL_SUPPORT_COLUMNS + [experimental_value_column]
    require_columns(df, required_columns, str(experimental_file))

    normalized_columns = list(metadata_column_map.values()) + EXPERIMENTAL_SUPPORT_COLUMNS
    df = normalize_string_columns(df, normalized_columns)
    df = df[required_columns].copy()
    df = expand_experimental_rows(df, experimental_value_column)
    df = add_log_experimental_column(df, experimental_value_column, experimental_log_value_column)
    check_unique_keys(df, DEFAULT_EXPERIMENTAL_ALIGN_COLUMNS, str(experimental_file))
    return df


def expand_experimental_rows(df: pd.DataFrame, experimental_value_column: str) -> pd.DataFrame:
    expanded_rows: list[dict[str, object]] = []
    for row in df.to_dict(orient="records"):
        target_chains = parse_csv_list(row.get("target_chain")) or parse_csv_list(row.get("chain_id")) or [""]
        keep_multiple = parse_csv_bool(row.get("keep_multiple_chains_in_struct"), default=False)
        for target_chain in target_chains:
            expanded_rows.append(
                {
                    "struct": get_struct_variant_name(
                        struct_name=row["struct"],
                        pdb_name=row["pdb_id"],
                        target_chain=target_chain,
                        keep_multiple_chains_in_struct=keep_multiple,
                    ),
                    "struct_base": row["struct"],
                    "target_chain": target_chain,
                    "ligname": row["ligname"],
                    "mutations": row["mutations"],
                    experimental_value_column: row[experimental_value_column],
                }
            )
    return pd.DataFrame(expanded_rows)


def add_log_experimental_column(
    df: pd.DataFrame,
    experimental_value_column: str,
    experimental_log_value_column: str,
) -> pd.DataFrame:
    df = df.copy()
    numeric_values = pd.to_numeric(df[experimental_value_column], errors="coerce")
    df[experimental_log_value_column] = np.where(numeric_values > 0, np.log10(numeric_values), np.nan)
    return df


def merge_calculation_with_experiment(
    calc_df: pd.DataFrame,
    experimental_df: pd.DataFrame,
    experimental_value_column: str,
) -> pd.DataFrame:
    merged_df = calc_df.merge(
        experimental_df,
        on=DEFAULT_EXPERIMENTAL_ALIGN_COLUMNS,
        how="inner",
        validate="many_to_one",
    )
    merged_df = merged_df.dropna(subset=[experimental_value_column]).reset_index(drop=True)
    return merged_df


def compute_correlation_metrics(
    df: pd.DataFrame,
    x_column: str,
    y_column: str,
) -> dict[str, float]:
    work_df = df[[x_column, y_column]].apply(pd.to_numeric, errors="coerce").dropna()
    n_points = len(work_df)

    if n_points < 2:
        return {
            "n_points": n_points,
            "pearson": math.nan,
            "spearman": math.nan,
        }

    x_values = work_df[x_column]
    y_values = work_df[y_column]
    pearson = x_values.corr(y_values, method="pearson")
    spearman = x_values.corr(y_values, method="spearman")
    return {
        "n_points": n_points,
        "pearson": pearson,
        "spearman": spearman,
    }


def build_segment_label(segment_values: tuple[object, ...], segment_columns: list[str]) -> str:
    return ", ".join(f"{column}={value}" for column, value in zip(segment_columns, segment_values))


def get_plot_df(
    df: pd.DataFrame,
    x_column: str,
    y_column: str,
    annotation_column: str = "mutations",
) -> pd.DataFrame:
    plot_columns = [x_column, y_column]
    if annotation_column in df.columns:
        plot_columns.append(annotation_column)
    plot_df = df[plot_columns].copy()
    plot_df[x_column] = pd.to_numeric(plot_df[x_column], errors="coerce")
    plot_df[y_column] = pd.to_numeric(plot_df[y_column], errors="coerce")
    return plot_df.dropna(subset=[x_column, y_column])


def plot_scatter(
    df: pd.DataFrame,
    x_column: str,
    y_column: str,
    title: str,
    output_path: Path,
    annotate_points: bool = False,
    annotation_column: str = "mutations",
) -> None:
    plot_df = get_plot_df(df, x_column, y_column, annotation_column=annotation_column)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(plot_df[x_column], plot_df[y_column], alpha=0.8)
    ax.set_xlabel(x_column)
    ax.set_ylabel(y_column)
    ax.set_title(title)
    ax.axhline(0.0, color="lightgray", linewidth=0.8)
    ax.axvline(0.0, color="lightgray", linewidth=0.8)
    ax.grid(alpha=0.25)

    if annotate_points and annotation_column in plot_df.columns:
        for _, row in plot_df.iterrows():
            ax.annotate(str(row[annotation_column]), (row[x_column], row[y_column]), fontsize=7, alpha=0.8)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def filter_segment_df(df: pd.DataFrame, segment_record: dict[str, object], segment_columns: list[str]) -> pd.DataFrame:
    segment_mask = pd.Series(True, index=df.index)
    for column in segment_columns:
        segment_mask &= df[column].eq(segment_record[column])
    return df.loc[segment_mask]


def analyze_segments(
    df: pd.DataFrame,
    calc_value_column: str,
    experimental_value_column: str,
    segment_columns: list[str],
    plot_prefix: str = "scatter",
) -> pd.DataFrame:
    summary_rows: list[dict[str, object]] = []

    for segment_values, segment_df in df.groupby(segment_columns, dropna=False):
        if not isinstance(segment_values, tuple):
            segment_values = (segment_values,)

        metrics = compute_correlation_metrics(segment_df, calc_value_column, experimental_value_column)
        summary_row = {
            "analysis_level": plot_prefix,
            **{column: value for column, value in zip(segment_columns, segment_values)},
            **metrics,
        }
        summary_rows.append(summary_row)

    if not summary_rows:
        return pd.DataFrame(
            columns=[
                "analysis_level",
                *segment_columns,
                "n_points",
                "pearson",
                "spearman",
                "abs_pearson",
                "abs_spearman",
            ]
        )

    summary_df = pd.DataFrame(summary_rows)
    summary_df["abs_pearson"] = summary_df["pearson"].abs()
    summary_df["abs_spearman"] = summary_df["spearman"].abs()
    return summary_df.sort_values(SUMMARY_SORT_COLUMNS, ascending=False).reset_index(drop=True)


def plot_segment_grid(
    df: pd.DataFrame,
    summary_df: pd.DataFrame,
    x_column: str,
    y_column: str,
    segment_columns: list[str],
    output_path: Path,
    subplots_per_row: int = 3,
    annotate_points: bool = False,
    annotation_column: str = "mutations",
) -> None:
    segment_records = summary_df[segment_columns + ["n_points", "pearson", "spearman"]].to_dict(orient="records")
    if not segment_records:
        return

    n_segments = len(segment_records)
    n_cols = max(1, subplots_per_row)
    n_rows = math.ceil(n_segments / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.8 * n_cols, 5.4 * n_rows), squeeze=False)
    axes_flat = axes.flatten()

    for ax, record in zip(axes_flat, segment_records):
        plot_df = get_plot_df(
            filter_segment_df(df, record, segment_columns),
            x_column,
            y_column,
            annotation_column=annotation_column,
        )
        segment_label = build_segment_label(tuple(record[column] for column in segment_columns), segment_columns)

        ax.scatter(plot_df[x_column], plot_df[y_column], alpha=0.8)
        ax.set_title(
            f"{segment_label}\npearson={record['pearson']:.3f}, spearman={record['spearman']:.3f}, n={record['n_points']}",
            fontsize=10,
        )
        ax.set_xlabel(x_column)
        ax.set_ylabel(y_column)
        ax.axhline(0.0, color="lightgray", linewidth=0.8)
        ax.axvline(0.0, color="lightgray", linewidth=0.8)
        ax.grid(alpha=0.25)

        if annotate_points and annotation_column in plot_df.columns:
            for _, row in plot_df.iterrows():
                ax.annotate(str(row[annotation_column]), (row[x_column], row[y_column]), fontsize=6, alpha=0.8)

    for ax in axes_flat[n_segments:]:
        ax.axis("off")

    fig.subplots_adjust(hspace=0.6, wspace=0.35, top=0.95, bottom=0.08)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_segment_individual(
    df: pd.DataFrame,
    summary_df: pd.DataFrame,
    x_column: str,
    y_column: str,
    segment_columns: list[str],
    output_dir: Path,
    output_prefix: str,
    plot_prefix: str,
    annotate_points: bool = False,
    annotation_column: str = "mutations",
) -> None:
    segment_records = summary_df[segment_columns + ["n_points", "pearson", "spearman"]].to_dict(orient="records")
    for record in segment_records:
        segment_df = filter_segment_df(df, record, segment_columns)
        filename_stub = "__".join(
            [plot_prefix] + [f"{column}-{str(record[column]).replace(' ', '_')}" for column in segment_columns]
        )
        title = (
            f"{build_segment_label(tuple(record[column] for column in segment_columns), segment_columns)}\n"
            f"pearson={record['pearson']:.3f}, spearman={record['spearman']:.3f}, n={record['n_points']}"
        )
        plot_scatter(
            df=segment_df,
            x_column=x_column,
            y_column=y_column,
            title=title,
            output_path=output_dir / f"{output_prefix}{filename_stub}.png",
            annotate_points=annotate_points,
            annotation_column=annotation_column,
        )


def plot_all_conditions_overlay(
    df: pd.DataFrame,
    summary_df: pd.DataFrame,
    x_column: str,
    y_column: str,
    segment_columns: list[str],
    output_path: Path,
) -> None:
    if summary_df.empty:
        return

    fig, ax = plt.subplots(figsize=(11, 7))
    cmap = plt.get_cmap("tab20")
    segment_records = summary_df[segment_columns].drop_duplicates().to_dict(orient="records")

    for idx, record in enumerate(segment_records):
        plot_df = get_plot_df(filter_segment_df(df, record, segment_columns), x_column, y_column)
        label = build_segment_label(tuple(record[column] for column in segment_columns), segment_columns)
        ax.scatter(plot_df[x_column], plot_df[y_column], alpha=0.7, color=cmap(idx % 20), label=label)

    ax.set_xlabel(x_column)
    ax.set_ylabel(y_column)
    ax.axhline(0.0, color="lightgray", linewidth=0.8)
    ax.axvline(0.0, color="lightgray", linewidth=0.8)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    fig.subplots_adjust(right=0.68)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def average_over_chains(
    merged_df: pd.DataFrame,
    calc_value_column: str,
    segment_columns: list[str],
    experimental_metric_column: str,
) -> pd.DataFrame:
    group_columns = [column for column in segment_columns if column != "target_chain"]
    group_columns += ["struct_base", "ligname", "mutations"]
    group_columns = list(dict.fromkeys(group_columns))

    averaged_df = (
        merged_df.groupby(group_columns, dropna=False, as_index=False)
        .agg(
            **{
                calc_value_column: (calc_value_column, "mean"),
                experimental_metric_column: (experimental_metric_column, "first"),
                "n_chains_averaged": (calc_value_column, "size"),
            }
        )
    )
    return averaged_df


def write_output_tables(
    output_dir: Path,
    output_tables: dict[str, pd.DataFrame],
    output_prefix: str = "",
    save_intermediate_tables: bool = False,
) -> None:
    output_filenames = dict(CORE_OUTPUT_FILENAMES)
    if save_intermediate_tables:
        output_filenames.update(INTERMEDIATE_OUTPUT_FILENAMES)
    for key, filename in output_filenames.items():
        output_tables[key].to_csv(output_dir / f"{output_prefix}{filename}", index=False)


def save_analysis_plots(
    df: pd.DataFrame,
    summary_df: pd.DataFrame,
    x_column: str,
    y_column: str,
    segment_columns: list[str],
    output_dir: Path,
    output_prefix: str,
    plot_prefix: str,
    plot_mode: str,
    subplots_per_row: int,
    annotate_points: bool,
) -> None:
    if plot_mode == "grid":
        plot_segment_grid(
            df=df,
            summary_df=summary_df,
            x_column=x_column,
            y_column=y_column,
            segment_columns=segment_columns,
            output_path=output_dir / f"{output_prefix}{plot_prefix}_scatter_grid.png",
            subplots_per_row=subplots_per_row,
            annotate_points=annotate_points,
        )
    else:
        plot_segment_individual(
            df=df,
            summary_df=summary_df,
            x_column=x_column,
            y_column=y_column,
            segment_columns=segment_columns,
            output_dir=output_dir,
            output_prefix=output_prefix,
            plot_prefix=plot_prefix,
            annotate_points=annotate_points,
        )

    plot_all_conditions_overlay(
        df=df,
        summary_df=summary_df,
        x_column=x_column,
        y_column=y_column,
        segment_columns=segment_columns,
        output_path=output_dir / f"{output_prefix}{plot_prefix}_scatter_overlay.png",
    )


def run_segment_analysis(
    df: pd.DataFrame,
    calc_value_column: str,
    experimental_value_column: str,
    segment_columns: list[str],
    plot_prefix: str,
) -> pd.DataFrame:
    return analyze_segments(
        df=df,
        calc_value_column=calc_value_column,
        experimental_value_column=experimental_value_column,
        segment_columns=segment_columns,
        plot_prefix=plot_prefix,
    )


def run_analysis(
    data_folder: str,
    data_subfolder: str,
    calc_value_columns: list[str],
    analysis_value_column: str,
    experimental_value_column: str,
    experimental_log_value_column: str,
    merge_id_columns: list[str],
    segment_columns: list[str],
    result_filenames: list[str],
    experimental_filename: str = "PA-NA_benchmark_wIC50foldchange.csv",
    output_subdir: str = "analysis",
    output_prefix: str = "",
    make_plots: bool = True,
    annotate_points: bool = False,
    save_intermediate_tables: bool = False,
    plot_mode: str = "grid",
    subplots_per_row: int = 3,
) -> dict[str, pd.DataFrame]:
    if plot_mode not in PLOT_MODES:
        raise ValueError(f"plot_mode must be one of {sorted(PLOT_MODES)}")

    base_dir = Path(data_folder)
    results_dir = base_dir / "yasara" / "Output" / data_subfolder
    experimental_dir = base_dir / "expdata" / data_subfolder
    output_dir = results_dir / output_subdir
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_plot_dir = output_dir / "raw_scatter"
    averaged_plot_dir = output_dir / "averaged_scatter"
    if make_plots:
        raw_plot_dir.mkdir(parents=True, exist_ok=True)
        averaged_plot_dir.mkdir(parents=True, exist_ok=True)

    result_files = [results_dir / filename for filename in result_filenames]
    missing_result_files = [str(file_path) for file_path in result_files if not file_path.exists()]
    if missing_result_files:
        raise FileNotFoundError(f"Missing result files: {missing_result_files}")

    if analysis_value_column not in calc_value_columns:
        raise ValueError(
            f"analysis_value_column={analysis_value_column!r} must be present in calc_value_columns={calc_value_columns!r}"
        )

    result_frames = [
        load_calculation_file(
            file_path=file_path,
            merge_id_columns=merge_id_columns,
            value_columns=calc_value_columns,
        )
        for file_path in result_files
    ]
    calc_long_df = concatenate_calculation_results_long(result_frames)
    calc_wide_df = merge_calculation_results_wide(result_frames, merge_id_columns, calc_value_columns)

    experimental_df = load_experimental_data(
        experimental_file=experimental_dir / experimental_filename,
        experimental_value_column=experimental_value_column,
        experimental_log_value_column=experimental_log_value_column,
    )
    merged_df = merge_calculation_with_experiment(
        calc_df=calc_long_df,
        experimental_df=experimental_df,
        experimental_value_column=experimental_value_column,
    )

    raw_summary_df = run_segment_analysis(
        df=merged_df,
        calc_value_column=analysis_value_column,
        experimental_value_column=experimental_log_value_column,
        segment_columns=segment_columns,
        plot_prefix="raw",
    )

    averaged_df = average_over_chains(
        merged_df=merged_df,
        calc_value_column=analysis_value_column,
        segment_columns=segment_columns,
        experimental_metric_column=experimental_log_value_column,
    )
    averaged_summary_df = run_segment_analysis(
        df=averaged_df,
        calc_value_column=analysis_value_column,
        experimental_value_column=experimental_log_value_column,
        segment_columns=segment_columns,
        plot_prefix="chain_averaged",
    )

    overall_summary_df = pd.concat([raw_summary_df, averaged_summary_df], axis=0, ignore_index=True, sort=False)
    overall_summary_df = overall_summary_df.sort_values(
        ["analysis_level", *SUMMARY_SORT_COLUMNS],
        ascending=[True, False, False],
    ).reset_index(drop=True)

    outputs = {
        "calc_wide": calc_wide_df,
        "calc_long": calc_long_df,
        "experimental": experimental_df,
        "merged_raw": merged_df,
        "merged_chain_averaged": averaged_df,
        "summary_raw": raw_summary_df,
        "summary_chain_averaged": averaged_summary_df,
        "summary_overall": overall_summary_df,
    }
    write_output_tables(
        output_dir,
        outputs,
        output_prefix=output_prefix,
        save_intermediate_tables=save_intermediate_tables,
    )
    if make_plots:
        save_analysis_plots(
            df=merged_df,
            summary_df=raw_summary_df,
            x_column=analysis_value_column,
            y_column=experimental_log_value_column,
            segment_columns=segment_columns,
            output_dir=raw_plot_dir,
            output_prefix=output_prefix,
            plot_prefix="raw",
            plot_mode=plot_mode,
            subplots_per_row=subplots_per_row,
            annotate_points=annotate_points,
        )
        save_analysis_plots(
            df=averaged_df,
            summary_df=averaged_summary_df,
            x_column=analysis_value_column,
            y_column=experimental_log_value_column,
            segment_columns=segment_columns,
            output_dir=averaged_plot_dir,
            output_prefix=output_prefix,
            plot_prefix="chain_averaged",
            plot_mode=plot_mode,
            subplots_per_row=subplots_per_row,
            annotate_points=annotate_points,
        )
    return outputs


if __name__ == "__main__":
    # Main analysis inputs: change these directly in your IDE when exploring.
    data_folder = "../influenza-resistance/"
    data_subfolder = ""

    # Result files to compare.
    result_filenames = [
        "PA-NA_benchmark_yasara2026_linux.csv",
        "PA-NA_benchmark_yasara2025_linux.csv",
    ]

    # Keep these columns when aligning files horizontally. ``target_chain`` is
    # included by default because ``setname`` does not necessarily distinguish
    # chains. If you need a chain-agnostic alignment for legacy files, remove
    # it here or switch to a different key such as ``['setname']``.
    merge_id_columns = DEFAULT_CALC_ALIGN_COLUMNS.copy()
    calc_value_columns = DEFAULT_CALC_VALUE_COLUMNS.copy()
    analysis_value_column = "ebindDDG"

    # Experimental data settings.
    experimental_filename = "PA-NA_benchmark_wIC50foldchange.csv"
    experimental_value_column = DEFAULT_EXPERIMENTAL_VALUE_COLUMN
    experimental_log_value_column = DEFAULT_EXPERIMENTAL_LOG_VALUE_COLUMN

    # Segmentation for plots and summary comparisons.
    segment_columns = DEFAULT_SEGMENT_COLUMNS.copy()

    # Output and plotting options.
    output_subdir = "analysis"
    output_prefix = "PA-NA_benchmark_"
    make_plots = True
    annotate_points = False
    save_intermediate_tables = False
    plot_mode = "grid" # "individual" #
    subplots_per_row = 3

    outputs = run_analysis(
        data_folder=data_folder,
        data_subfolder=data_subfolder,
        calc_value_columns=calc_value_columns,
        analysis_value_column=analysis_value_column,
        experimental_value_column=experimental_value_column,
        experimental_log_value_column=experimental_log_value_column,
        merge_id_columns=merge_id_columns,
        segment_columns=segment_columns,
        result_filenames=result_filenames,
        experimental_filename=experimental_filename,
        output_subdir=output_subdir,
        output_prefix=output_prefix,
        make_plots=make_plots,
        annotate_points=annotate_points,
        save_intermediate_tables=save_intermediate_tables,
        plot_mode=plot_mode,
        subplots_per_row=subplots_per_row,
    )

    print("Saved analysis outputs to:", Path(data_folder) / "yasara" / "Output" / data_subfolder / output_subdir)
    print(outputs["summary_overall"].to_string(index=False))
