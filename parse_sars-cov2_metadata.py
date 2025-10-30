"""Command-line utility to reproduce the processing done in tercera_ola.ipynb.

The original notebook orchestrated the SARS-CoV-2 data harmonisation for the
Andalusian sequencing effort.  This module keeps the same steps in a more
maintainable Python script so it can be executed without Jupyter.

Usage::

    python tercera_ola.py

The script assumes the same directory layout that the notebook expected.  You
can override the base directory by exporting ``COVID19_SSPA_BASE``.
"""
from __future__ import annotations

import datetime as dt
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd


BASE_PATH = Path(os.environ.get("COVID19_SSPA_BASE", "/mnt/NAS/projects/virus_SSPA/COVID19_SSPA/"))
NEXTSTRAIN_DATA = BASE_PATH / "nextstrain_all/data"
WAVE1_PATH = BASE_PATH / "other_data/1wave_data"
SSP_PATH = BASE_PATH / "secuenciacion_salud_publica"
METADATA_PATH = SSP_PATH / "metadata"
DATA_DIR = Path("./data")
CONFIG_DIR = Path("./config")
WEB_REPORTS_DIR = SSP_PATH / "analysis_MN908947/reports/web"


def fasta_get_by_name(fasta: Path | str, ids: Iterable[str], delimit: str = " ", column: int = 0) -> Dict[str, str]:
    """Return a ``header -> sequence`` dictionary for the requested identifiers."""

    fasta_path = Path(fasta)
    targets = set(ids)
    if not fasta_path.exists():
        raise FileNotFoundError(f"Missing FASTA file: {fasta_path}")

    sequences: Dict[str, str] = {}
    current_header = ""
    with fasta_path.open("r") as handle:
        for raw_line in handle:
            line = raw_line.rstrip()
            if line.startswith(">"):
                current_header = line[1:]
                sequences[current_header] = ""
                continue
            sequences[current_header] += line

    subset: Dict[str, str] = {}
    for header, sequence in sequences.items():
        token = header.split(delimit)[column]
        if token in targets:
            subset[header] = sequence
    return subset


def read_batch_plan() -> pd.DataFrame:
    plan_file = SSP_PATH / "sarscov2_sequencing_analysis_info.csv"
    return pd.read_csv(plan_file, low_memory=False)


def clean_lab_ids(series: pd.Series) -> pd.Series:
    cleaned = series.copy()
    cleaned = cleaned.str.replace(r"^1442", "", regex=True)
    cleaned = cleaned.str.replace(r"^1012", "", regex=True)
    cleaned = cleaned.str.replace(r"^47", "", regex=True)
    return cleaned


def build_samples(plan: pd.DataFrame) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    for run_id, qc_path, viralrecon_version in plan[["run_id", "path_QC", "viralrecon_version"]].values:
        qc = pd.read_csv(qc_path, dtype={"laboratoryID": str, "batch": str})
        qc.drop(columns=["qc_nextclade", "totalGaps_nextclade", "totalMutations_nextclade"], inplace=True, errors="ignore")
        qc["tecnología"] = "illumina" if viralrecon_version != "iontorrent" else viralrecon_version
        if run_id == "AND_2_husc_010221":
            qc["laboratoryID"] = qc["laboratoryID"].apply(lambda value: re.sub(r"^1442", "", value))
        frames.append(qc)

    samples = pd.concat(frames, ignore_index=True)
    samples["laboratoryID"] = clean_lab_ids(samples["laboratoryID"].astype(str))
    return samples


def build_mapping(plan: pd.DataFrame) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    for run_id in plan["run_id"].values:
        location = ""
        if "huvr" in run_id:
            location = "huvr"
        elif "husc" in run_id:
            location = "husc"
        elif "huvn" in run_id:
            location = "huvn"
        elif "hrum" in run_id:
            location = "hrum"
        else:
            continue

        mapping_file = METADATA_PATH / location / run_id / f"{run_id}_sampleIDlabo_mapping_andID.csv"
        lab = pd.read_csv(mapping_file, header=None, dtype={2: str})
        if run_id == "AND_2_husc_010221":
            lab[2] = lab[2].apply(lambda value: re.sub(r"^1442", "", value))
        frames.append(lab)

    mapping = pd.concat(frames, ignore_index=True)
    mapping.columns = ["sample", "and_id", "num_lab"]
    mapping = mapping.loc[~mapping.duplicated(keep="first")]  # keep first duplicate to avoid confusion
    mapping.reset_index(drop=True, inplace=True)
    mapping["num_lab"] = clean_lab_ids(mapping["num_lab"].astype(str))
    return mapping


def build_lineages(plan: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for run_id, analysis_path in plan[["run_id", "analysis_path"]].values:
        pangolin_file = Path(analysis_path) / "pangolin" / f"{run_id}_pangolin.csv"
        frames.append(pd.read_csv(pangolin_file))

    wave1_file = WAVE1_PATH / "pangolin/andalusia_seq_pangolin.csv"
    frames.append(pd.read_csv(wave1_file))
    return pd.concat(frames, ignore_index=True)


def build_nextclade(plan: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for analysis_path, run_id in plan[["analysis_path", "run_id"]].values:
        csv_file = Path(analysis_path) / "nextclade" / f"{run_id}.csv"
        frames.append(pd.read_csv(csv_file, sep=";"))

    wave1_file = WAVE1_PATH / "nextclade/andalusia_seq_nextclade.csv"
    frames.append(pd.read_csv(wave1_file, sep=";"))
    nextclade = pd.concat(frames, ignore_index=True)
    return nextclade[["seqName", "clade", "clade_display", "clade_who", "Nextclade_pango"]]


def build_metadata() -> pd.DataFrame:
    dtype = {
        "IDlaboratorio": str,
        "Código Postal": str,
        "Edad (años)": "Int64",
        "Nº Contactos identificados": str,
        "Identificador": "Int64",
        "Tipo": str,
    }
    remove_cols = [
        "Evolución",
        "Fecha Alta (dd/mm/aaaa)",
        "Fecha de defunción",
        "Fue Hospitalizado",
        "Hospital de ingreso",
        "Identificador",
        "Identificador Brote asociado",
        "Tipo",
        "Ámbito",
        "Clasificación reinfección",
        "Otros cuadros respiratorios graves",
        "Asintomático (sin antecedentes de síntomas)",
        "Categoría profesional",
        "Contacto estrecho",
        "Diabetes",
        "Embarazo",
        "Enfermedad cardiovascular",
        "Enfermedad pulmonar crónica",
        "Estancia en UCI",
        "Factores de Riesgo y Enfermedad de Base COVID-19",
        "Fallo renal agudo",
        "Fecha alta UCI",
        "Fecha ingreso UCI",
        "Hipertensión arterial (HTA)",
        "Inmunosupresión",
        "Institucionalizados (revisar ámbito)",
        "Neumonía",
        "Otros factores",
        "Otros factores, especificar",
        "Síndrome Distress respiratorio",
        "Ventilación mecánica",
        "Cáncer",
        "Estudio Contactos",
        "Fecha primera vacuna",
        "Fecha segunda vacuna",
        "Nº Contactos confirmados como caso",
        "Nº Contactos identificados",
        "Nombre primera vacuna",
        "Nombre segunda vacuna",
        "FERE1",
        "RE1",
        "FERE2",
        "RE2",
        "FERE3",
        "RE3",
        "DATABASE RLT2_NUHSA",
    ]

    originals = [
        pd.read_csv(METADATA_PATH / "huvr/AND_1_huvr_050121/AND_1_huvr_050121_sp.csv", sep="|", dtype=dtype),
        pd.read_csv(METADATA_PATH / "husc/AND_2_husc_010221/AND_2_husc_010221_sp.csv", sep="|", dtype=dtype),
        pd.read_csv(METADATA_PATH / "huvr/AND_3_huvr_020221/AND_3_huvr_020221_sp.csv", sep="|", dtype=dtype),
        pd.read_csv(METADATA_PATH / "husc/AND_4_husc_200221/AND_4_husc_200221_sp.csv", sep="|", dtype=dtype),
        pd.read_csv(METADATA_PATH / "huvr/AND_5_huvr_240221/AND_5_huvr_240221_sp.csv", sep="|", dtype=dtype),
        pd.read_csv(METADATA_PATH / "huvr/AND_6_huvr_270221/AND_6_huvr_270221_sp.csv", sep="|", dtype=dtype),
        pd.read_csv(METADATA_PATH / "huvr/AND_7_huvr_080321/AND_7_huvr_080321_sp.csv", sep=",", dtype=dtype),
    ]
    for frame in originals:
        frame["ori"] = "ori_met"

    latest = [
        pd.read_csv(METADATA_PATH / "huvr/huvr_acumulativo_latest_sp.csv", sep="|", low_memory=False, dtype=dtype).assign(ori="met_huvr"),
        pd.read_csv(METADATA_PATH / "husc/husc_acumulativo_latest_sp.csv", sep="|", low_memory=False, dtype=dtype).assign(ori="met_husc"),
        pd.read_csv(METADATA_PATH / "huvn/huvn_acumulativo_latest_sp.csv", sep="|", low_memory=False, dtype=dtype).assign(ori="met_huvn"),
        pd.read_csv(METADATA_PATH / "hrum/hrum_acumulativo_latest.csv", sep="|", low_memory=False, dtype=dtype).assign(ori="met_hrum"),
    ]

    original_df = pd.concat(originals, ignore_index=True)
    latest_df = pd.concat(latest, ignore_index=True)
    original_df = original_df.loc[
        (~original_df["IDlaboratorio"].isin(latest_df["IDlaboratorio"]))
        & (~original_df["NUHSA"].isin(latest_df["NUHSA"]))
    ]

    metadata = pd.concat([latest_df, original_df], ignore_index=True)
    metadata.drop(columns=remove_cols, inplace=True)
    metadata.drop_duplicates(inplace=True, keep="first")
    metadata.reset_index(inplace=True, drop=False)

    duplicated_ids = metadata.loc[metadata.duplicated(["IDlaboratorio"], keep=False)]
    if not duplicated_ids.empty:
        print(f"Warning: {duplicated_ids.shape[0]} duplicated laboratory IDs in metadata")
        print(duplicated_ids[["IDlaboratorio", "NUHSA"]].sort_values("IDlaboratorio").head())
        duplicated_ids.sort_values("IDlaboratorio").to_excel(
            "metadata_laboratory_id_duplicated.xlsx", index=False
        )
    metadata.drop_duplicates("IDlaboratorio", inplace=True, keep="first")

    needs_padding = metadata["Código Postal"].str.len() == 4
    metadata.loc[needs_padding & metadata["Código Postal"].str.match(r"^4"), "Código Postal"] = (
        "0" + metadata.loc[needs_padding & metadata["Código Postal"].str.match(r"^4"), "Código Postal"]
    )
    metadata["Municipio"] = metadata["Municipio"].str.strip()
    return metadata


def merge_metadata(samples: pd.DataFrame, lineage: pd.DataFrame, nextclade: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    merged = samples.copy()
    print(f"Total samples: {merged.shape[0]}")
    merged = merged.merge(lineage, how="left", left_on="sample", right_on="taxon")
    print(f"Samples merged with lineage: {merged.shape[0]}")
    merged = merged.merge(nextclade, how="left", left_on="sample", right_on="seqName")
    print(f"Samples merged with nextclade: {merged.shape[0]}")
    merged = merged.merge(metadata, how="left", left_on="laboratoryID", right_on="IDlaboratorio")
    print(f"Samples merged with metadata: {merged.shape[0]}")
    merged.drop_duplicates(inplace=True, keep="first")
    print(f"Samples without duplicates: {merged.shape[0]}")
    duplicated_lab = merged.loc[merged["laboratoryID"].duplicated(keep=False)]
    if not duplicated_lab.empty:
        print(f"Warning: {duplicated_lab.shape[0]} duplicated laboratory IDs in samples")
        duplicated_lab.sort_values("laboratoryID").to_excel("samples_laboratory_id_duplicated.xlsx", index=False)
    return merged


def assign_geography(samples: pd.DataFrame) -> pd.DataFrame:
    enriched = samples.copy()
    enriched["origen"] = "SAS"
    enriched.loc[enriched["Código Postal"].isnull(), "Código Postal"] = "0"

    postal_to_province = {
        "41": "Sevilla",
        "14": "Córdoba",
        "21": "Huelva",
        "23": "Jaén",
        "18": "Granada",
        "29": "Málaga",
        "04": "Almería",
        "11": "Cádiz",
    }
    for prefix, province in postal_to_province.items():
        enriched.loc[enriched["Código Postal"].str.startswith(prefix), "provincia"] = province

    enriched["Municipio"] = enriched["Municipio"].str.strip()

    location_dict = DATA_DIR / "location_names_dict.tsv"
    if location_dict.exists():
        with location_dict.open("r") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                src, dst = line.split("\t")
                enriched.loc[enriched["Municipio"] == src, "Municipio"] = dst

    lat_long_file = CONFIG_DIR / "lat_long.tsv"
    if lat_long_file.exists():
        lat_long = pd.read_csv(lat_long_file, sep="\t", names=["location", "name", "lat", "long"])
        missing = enriched.loc[
            (~enriched["Municipio"].isin(lat_long["name"])) & enriched["Municipio"].notnull(),
            "Municipio",
        ]
        print("Municipalities without coordinates")
        print(missing)

    province_dict = DATA_DIR / "location_province_dict.tsv"
    if province_dict.exists():
        with province_dict.open("r") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                location, province = line.split("\t")
                enriched.loc[enriched["Municipio"] == location, "provincia"] = province

    print(
        "Municipalities without province\n",
        enriched.loc[(enriched["provincia"].isnull()) & enriched["Municipio"].notnull(), "Municipio"].value_counts(),
    )
    return enriched


def assign_hospitals(samples: pd.DataFrame) -> pd.DataFrame:
    enriched = samples.copy()
    enriched.loc[enriched["batch"].str.contains("huvr"), "hospital de referencia"] = "HUVR"
    enriched.loc[enriched["batch"].str.contains("husc"), "hospital de referencia"] = "HUSC"

    province_to_hospital = {
        "Sevilla": "Hospital Universitario Virgen del Rocío",
        "Córdoba": "Hospital Universitario Reina Sofía",
        "Huelva": "Hospital Universitario Juan Ramón Jiménez",
        "Jaén": "Hospital Universitario de Jaén",
        "Granada": "Hospital Universitario San Cecilio",
        "Málaga": "Hospital Universitario Virgen de la Victoria",
        "Almería": "Hospital Universitario Torrecárdenas",
        "Cádiz": "Hospital Universitario Puerta del Mar",
    }
    for province, hospital in province_to_hospital.items():
        enriched.loc[enriched["provincia"] == province, "hospital"] = hospital

    enriched.loc[enriched["batch"].str.contains("huvn"), "hospital"] = "Hospital Universitario Vigen de las Nieves"
    enriched.loc[enriched["batch"].str.contains("hurm"), "hospital"] = "Hospital Universitario Regional de Málaga"
    enriched.drop_duplicates(keep="first", inplace=True)
    enriched.loc[enriched["seqPlatform"].isnull(), "seqPlatform"] = "Illumina"
    enriched.loc[
        enriched["batch"].str.contains("huvr"), "centro_secuenciacion"
    ] = "microbiología HU Virgen del Rocio"
    return enriched


def update_dates(samples: pd.DataFrame) -> pd.DataFrame:
    result = samples.copy()
    result["fecha_nextstrain"] = result["Fecha Inicio Síntomas (dd/mm/aaaa)"]
    result.loc[result["fecha_nextstrain"].isnull(), "fecha_nextstrain"] = result["Fecha del caso"]
    result.loc[result["fecha_nextstrain"].isnull(), "fecha_nextstrain"] = result["Fecha_diag"]
    result["fecha_nextstrain"] = pd.to_datetime(result["fecha_nextstrain"], format="%d/%m/%Y")

    occident_path = METADATA_PATH / "huvr/MUESTRAS_SECUENCIACION_SEMANALES_acumulada_latest.csv"
    occident = pd.read_csv(occident_path, sep=";")
    occident_errors = occident[["Identificacion", "NUHSA", "Fecha"]].loc[
        pd.to_datetime(occident["Fecha"], format="%d/%m/%Y", errors="coerce").isnull()
    ]
    if occident_errors.shape[0] > 0:
        print(">>> Dates with wrong format (HUVR)")
        print(occident_errors)
        occident["Fecha"] = pd.to_datetime(occident["Fecha"], format="%d/%m/%Y", errors="coerce")
    else:
        print(">>> HUVR dates format OK")
        occident["Fecha"] = pd.to_datetime(occident["Fecha"], format="%d/%m/%Y")

    today = dt.datetime.now()
    future = occident.loc[occident["Fecha"] > today]
    if not future.empty:
        print(">>> Future dates in HUVR")
        print(future[["Identificacion", "NUHSA", "Fecha"]])
        occident_errors = pd.concat([occident_errors, future], ignore_index=True)
        occident = occident.loc[occident["Fecha"] <= today]
    else:
        print(">>> HUVR future dates OK")

    old = occident.loc[occident["Fecha"] < "2020-01-01"]
    if not old.empty:
        print(">>> Old dates in HUVR")
        print(old[["Identificacion", "NUHSA", "Fecha"]])
        occident_errors = pd.concat([occident_errors, old], ignore_index=True)
        occident = occident.loc[occident["Fecha"] >= "2020-01-01"]
    else:
        print(">>> HUVR old dates OK")

    occident_errors.to_csv("date_occi_error.tsv", sep="\t", index=False)
    occident = occident.drop_duplicates(subset=["Fecha", "Identificacion"], keep="first")
    occident = occident.drop_duplicates(subset=["NUHSA", "Identificacion"], keep="first")

    result = result.merge(occident[["Fecha", "Identificacion"]], how="left", left_on="laboratoryID", right_on="Identificacion")
    result.rename(columns={"Fecha": "Fecha_micro_huvr"}, inplace=True)
    result["Fecha_micro"] = pd.to_datetime(result["Fecha_micro_huvr"], format="%Y-%m-%d", errors="coerce")

    prison_ids: List[str] = []
    if "Justificacion clinica" in occident.columns:
        prison_cases = occident.loc[occident["Justificacion clinica"] == "Brote prisión", ["Identificacion"]]
        prison_ids = prison_cases["Identificacion"].tolist()
    result.loc[result["laboratoryID"].isin(prison_ids), "Municipio"] = "Córdoba"
    result.loc[result["laboratoryID"].isin(prison_ids), "provincia"] = "Córdoba"

    orient_path = SSP_PATH / "metadata/husc/Actualizacion_Salud_Publica_latest.csv"
    orient = pd.read_csv(orient_path, sep=";", usecols=["Identificador", "NUHSA", "Fecha de la pueba"])
    orient.rename(columns={"Identificador": "ori_temp_id"}, inplace=True)
    orient_errors = orient.loc[orient["Fecha de la pueba"].isnull(), ["ori_temp_id", "NUHSA", "Fecha de la pueba"]]
    orient = orient.loc[orient["Fecha de la pueba"].notnull()]
    orient = orient.loc[orient["Fecha de la pueba"].str.contains("/")]
    orient["ori_temp_id"] = clean_lab_ids(orient["ori_temp_id"].astype(str))

    for suffix in ("2020", "2021", "2022", "2023", "2024", "2025"):
        orient["Fecha de la pueba"] = orient["Fecha de la pueba"].str.replace(f"/{suffix}", f"/{suffix[-2:]}")

    orient[["dd", "mm", "yy"]] = orient["Fecha de la pueba"].str.split("/", expand=True)
    orient["mmint"] = orient["mm"].astype(int)
    orient.loc[orient["mmint"] > 12, "Fecha de la pueba"] = (
        orient["mm"].astype(str) + "-" + orient["dd"].astype(str) + "-" + orient["yy"].astype(str)
    )

    orient["ori_temp_date"] = orient["Fecha de la pueba"]
    parsed = pd.to_datetime(orient["ori_temp_date"], format="%d/%m/%y", errors="coerce")
    failures = orient.loc[parsed.isnull(), ["ori_temp_id", "ori_temp_date", "Fecha de la pueba"]]
    if not failures.empty:
        print(">>> HUSC dates with wrong format")
        print(failures)
        orient_errors = pd.concat([orient_errors, failures], ignore_index=True)
    orient["ori_temp_date"] = parsed

    future = orient.loc[orient["ori_temp_date"] > today]
    if not future.empty:
        print(">>> Future dates in HUSC")
        print(future[["ori_temp_id", "NUHSA", "Fecha de la pueba", "ori_temp_date"]])
        orient_errors = pd.concat([orient_errors, future], ignore_index=True)
        orient = orient.loc[orient["ori_temp_date"] <= today]
    else:
        print(">>> HUSC future dates OK")

    old = orient.loc[orient["ori_temp_date"] < "2020-01-01"]
    if not old.empty:
        print(">>> Old dates in HUSC")
        print(old[["ori_temp_id", "NUHSA", "Fecha de la pueba", "ori_temp_date"]])
        orient_errors = pd.concat([orient_errors, old], ignore_index=True)
        orient = orient.loc[orient["ori_temp_date"] >= "2020-01-01"]
    else:
        print(">>> HUSC old dates OK")

    orient_errors.to_csv("date_ori_errors.tsv", sep="\t", index=False)
    orient.drop_duplicates(subset=["ori_temp_date", "ori_temp_id"], keep="first", inplace=True)
    orient.drop_duplicates(subset=["NUHSA", "ori_temp_id"], keep="first", inplace=True)

    duplicates = orient.loc[orient["ori_temp_id"].duplicated(keep=False), "ori_temp_id"].unique()
    for lab_id in duplicates:
        if orient.loc[orient["ori_temp_id"] == lab_id, "NUHSA"].nunique() > 1:
            orient.drop(orient.loc[orient["ori_temp_id"] == lab_id].index, inplace=True)

    result = result.merge(orient[["ori_temp_date", "ori_temp_id"]], how="left", left_on="laboratoryID", right_on="ori_temp_id")
    result["fecha_micro_husc"] = result["ori_temp_date"]
    result.loc[result["Fecha_micro"].isnull() & result["ori_temp_date"].notnull(), "Fecha_micro"] = result["ori_temp_date"]

    malaga_path = SSP_PATH / "metadata/hrum/HRUM_INFORME_latest.csv"
    malaga = pd.read_csv(malaga_path, sep=";", dtype=str)[["Petición", "Fecha registro"]]
    malaga_errors = [malaga.loc[malaga["Fecha registro"].isnull()]]
    malaga_parsed = pd.to_datetime(malaga["Fecha registro"], format="%m/%d/%Y", errors="coerce")
    malaga_errors.append(malaga.loc[malaga_parsed.isnull()])
    malaga["Fecha registro"] = malaga_parsed
    future = malaga.loc[malaga["Fecha registro"] > today]
    if not future.empty:
        print(">>> Future dates in HRUM")
        print(future)
        malaga_errors.append(future)
        malaga = malaga.loc[malaga["Fecha registro"] <= today]
    else:
        print(">>> HRUM future dates OK")
    old = malaga.loc[malaga["Fecha registro"] < "2020-01-01"]
    if not old.empty:
        print(">>> Old dates in HRUM")
        print(old)
        malaga_errors.append(old)
        malaga = malaga.loc[malaga["Fecha registro"] >= "2020-01-01"]
    else:
        print(">>> HRUM old dates OK")

    result = result.merge(malaga, how="left", left_on="laboratoryID", right_on="Petición")
    result["fecha_micro_hurm"] = result["Fecha registro"]
    result.loc[result["Fecha_micro"].isnull() & result["Fecha registro"].notnull(), "Fecha_micro"] = result["Fecha registro"]

    granada_path = METADATA_PATH / "huvn/INFORME_ACUMULADO_SECUENCIACION_SARS-COV-2_actualizado_latest.csv"
    granada = pd.read_csv(granada_path, sep=";", dtype=str)
    granada["date_huvn"] = granada["Fecha de toma de muestra"].fillna(granada["Fecha de recepción de la muestra"])
    format_errors = granada.loc[pd.to_datetime(granada["date_huvn"], format="%d/%m/%Y", errors="coerce").isnull()]
    print(f">>> Errors in HUVN date format %d/%m/%Y: {format_errors.shape[0]}")
    granada["date_huvn"] = pd.to_datetime(granada["date_huvn"], format="%d/%m/%Y", errors="coerce")
    future = granada.loc[granada["date_huvn"] > today]
    if not future.empty:
        print(">>> Future dates in HUVN")
        print(future)
        granada = granada.loc[granada["date_huvn"] <= today]
    else:
        print(">>> HUVN future dates OK")
    old = granada.loc[granada["date_huvn"] < "2020-01-01"]
    if not old.empty:
        print(">>> Old dates in HUVN")
        print(old)
        granada = granada.loc[granada["date_huvn"] >= "2020-01-01"]
    else:
        print(">>> HUVN old dates OK")

    result = result.merge(granada[["date_huvn", "Número"]], how="left", left_on="laboratoryID", right_on="Número")
    result["fecha_micro_huvn"] = result["date_huvn"]
    result.loc[result["Fecha_micro"].isnull() & result["date_huvn"].notnull(), "Fecha_micro"] = result["date_huvn"]

    result.loc[result["Fecha_micro"].notnull(), "fecha_nextstrain"] = result["Fecha_micro"]
    result.drop(columns=["ori_temp_date", "ori_temp_id", "Petición", "Fecha registro", "Número", "date_huvn"], inplace=True)

    return result


def build_auspice(samples: pd.DataFrame, nextclade: pd.DataFrame, lineage: pd.DataFrame) -> pd.DataFrame:
    subset = samples[
        [
            "sample",
            "laboratoryID",
            "NUHSA",
            "fecha_nextstrain",
            "País del caso",
            "provincia",
            "Municipio",
            "origen",
            "hospital",
            "hospital de referencia",
            "Centro de salud",
            "centro_secuenciacion",
            "lineage",
            "clade",
            "clade_display",
            "nextstrain_genome",
            "selected_for_nextstrain",
            "batch",
            "percentage_reference_genome_in_consensus_ivar",
            "seqPlatform",
        ]
    ].copy()

    subset.rename(columns={"sample": "strain", "fecha_nextstrain": "date", "Municipio": "location"}, inplace=True)
    subset["source"] = "SAS"
    subset["location"] = subset["location"].str.strip()

    location_names = DATA_DIR / "location_names_dict.tsv"
    if location_names.exists():
        translation: Dict[str, str] = {}
        with location_names.open("r") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                origin, target = line.split("\t")
                translation[origin] = target
        for origin, target in translation.items():
            subset.loc[:, "location"] = subset["location"].replace(origin, target)

    print(f"Auspice records: {subset.shape[0]}")
    selected = subset.loc[subset["selected_for_nextstrain"] == "yes"].copy()
    print(f"Auspice selected for nextstrain: {selected.shape[0]}")

    first_wave = pd.read_csv(WAVE1_PATH / "metadata.tsv", sep="\t")
    first_wave["1wave"] = "yes"
    first_wave["nextstrain_genome"] = np.nan
    first_wave["date"] = pd.to_datetime(first_wave["date"], format="%Y-%m-%d")
    first_wave.drop(columns=["clade", "lineage"], inplace=True)
    first_wave = first_wave.merge(nextclade[["seqName", "clade", "clade_display"]], how="left", left_on="strain", right_on="seqName")
    first_wave = first_wave.merge(lineage[["taxon", "lineage"]], how="left", left_on="strain", right_on="taxon")
    first_wave.drop(columns=["seqName", "taxon"], inplace=True)

    selected = pd.concat([selected, first_wave], ignore_index=True)
    print(f"Auspice with first wave: {selected.shape[0]}")

    selected = selected.loc[selected["date"].notnull()]
    print(f"Auspice with date: {selected.shape[0]}")

    alpha_err = selected.loc[(selected["date"] < "2020-12-01") & (selected["lineage"] == "B.1.1.7")]
    delta_err = selected.loc[(selected["date"] < "2021-05-01") & (selected["lineage"].str.startswith("AY"))]
    delta_err2 = selected.loc[(selected["date"] < "2021-05-01") & (selected["lineage"] == "B.1.617.2")]
    omicron_err = selected.loc[(selected["date"] < "2021-12-01") & (selected["lineage"].str.startswith("BA"))]
    err = pd.concat([alpha_err, delta_err, delta_err2, omicron_err])
    err.to_excel("error_date_lineage.xlsx", index=False)

    selected = selected.loc[~selected["strain"].isin(err["strain"])]
    print(f"Auspice after date-lineage filter: {selected.shape[0]}")

    outliers_file = NEXTSTRAIN_DATA / "remove_SARSCOV2_samples_outliers.txt"
    if outliers_file.exists():
        with outliers_file.open("r") as handle:
            outliers = [line.strip() for line in handle if line.strip()]
        selected = selected.loc[~selected["strain"].isin(outliers)]
    print(f"Auspice after AND id filter: {selected.shape[0]}")

    selected.drop_duplicates("strain", keep="first", inplace=True)
    print(f"Auspice without strain duplicates: {selected.shape[0]}")

    today = dt.datetime.now()
    selected = selected.loc[selected["date"] < today]
    print(f"Auspice without future dates: {selected.shape[0]}")

    selected.loc[selected["batch"].isnull(), "batch"] = ""

    WEB_REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = dt.datetime.now().strftime("%d%m%y")
    aggregated = selected[["strain", "batch", "date", "lineage", "location", "provincia", "hospital", "hospital de referencia"]].copy()
    aggregated.to_csv(WEB_REPORTS_DIR / f"agregado_web_noloc_{timestamp}.tsv", sep="\t", index=False)
    aggregated.loc[aggregated["batch"].str.contains("husc")].to_csv(
        WEB_REPORTS_DIR / f"agregado_web_noloc_husc_{timestamp}.tsv", sep="\t", index=False
    )
    aggregated.loc[aggregated["batch"].str.contains("huvn")].to_csv(
        WEB_REPORTS_DIR / f"agregado_web_noloc_huvn_{timestamp}.tsv", sep="\t", index=False
    )
    aggregated.loc[aggregated["batch"].str.contains("huvr")].to_csv(
        WEB_REPORTS_DIR / f"agregado_web_noloc_huvr_{timestamp}.tsv", sep="\t", index=False
    )

    selected = selected.loc[selected["location"].notnull()]
    print(f"Auspice with location: {selected.shape[0]}")
    selected = selected.loc[selected.get("source", "") != "SAS-imputed"]
    print(f"Auspice without SAS-imputed: {selected.shape[0]}")

    aggregated = selected[["strain", "batch", "date", "lineage", "location", "provincia", "hospital", "hospital de referencia"]]
    aggregated.to_csv(WEB_REPORTS_DIR / f"agregado_web_{timestamp}.tsv", sep="\t", index=False)
    aggregated.loc[aggregated["batch"].str.contains("husc")].to_csv(
        WEB_REPORTS_DIR / f"agregado_web_husc_{timestamp}.tsv", sep="\t", index=False
    )
    aggregated.loc[aggregated["batch"].str.contains("huvn")].to_csv(
        WEB_REPORTS_DIR / f"agregado_web_huvn_{timestamp}.tsv", sep="\t", index=False
    )
    aggregated.loc[aggregated["batch"].str.contains("huvr")].to_csv(
        WEB_REPORTS_DIR / f"agregado_web_huvr_{timestamp}.tsv", sep="\t", index=False
    )

    selected["date"] = pd.to_datetime(selected["date"]).dt.strftime("%Y-%m-%d")

    lat_long = pd.read_csv(CONFIG_DIR / "lat_long.tsv", sep="\t", header=None)
    lat_long_names = set(lat_long[1].values)
    locs = set(selected["location"].values)
    missing = locs - lat_long_names
    if missing:
        print(f"Missing coordinates for: {sorted(missing)}")

    four_months_ago = dt.date.today() - dt.timedelta(days=int(4 * 365 / 12))
    cutoff = pd.to_datetime(four_months_ago.strftime("%Y-%m-%d"))
    recent_mask = pd.to_datetime(selected["date"]) >= cutoff
    selected.loc[recent_mask, "include"] = "yes"

    aggregated_no_filter = samples[[
        "sample",
        "lineage",
        "laboratoryID",
        "batch",
        "hospital",
        "hospital de referencia",
        "fecha_nextstrain",
    ]].copy()
    aggregated_no_filter.to_csv(WEB_REPORTS_DIR / f"agregado_without_QC_{timestamp}.tsv", sep="\t", index=False)
    aggregated_no_filter.loc[aggregated_no_filter["batch"].str.contains("husc")].to_csv(
        WEB_REPORTS_DIR / f"agregado_without_QC_husc_{timestamp}.tsv", sep="\t", index=False
    )
    aggregated_no_filter.loc[aggregated_no_filter["batch"].str.contains("huvn")].to_csv(
        WEB_REPORTS_DIR / f"agregado_without_QC_huvn_{timestamp}.tsv", sep="\t", index=False
    )
    aggregated_no_filter.loc[aggregated_no_filter["batch"].str.contains("huvr")].to_csv(
        WEB_REPORTS_DIR / f"agregado_without_QC_huvr_{timestamp}.tsv", sep="\t", index=False
    )

    return selected


def write_sequences(auspice: pd.DataFrame, samples: pd.DataFrame) -> None:
    consensus_seq: Dict[str, str] = {}
    for strain, path in auspice[["strain", "nextstrain_genome"]].loc[auspice["nextstrain_genome"].notnull()].values:
        if path in {"unknown", "-"}:
            continue
        fasta_path = Path(path.replace("/analysis_MN90894/", "/analysis_MN908947/"))
        with fasta_path.open("r") as handle:
            key = ""
            for raw_line in handle:
                line = raw_line.strip()
                if line.startswith(">"):
                    key = strain
                    consensus_seq[key] = ""
                    continue
                consensus_seq[key] += line

    wave1_fasta = WAVE1_PATH / "andalusia_seq.fasta"
    ids_wave1 = auspice.loc[auspice.get("1wave", "") == "yes", "strain"].values
    wave1_seq = fasta_get_by_name(wave1_fasta, ids_wave1)

    seq = {**consensus_seq, **wave1_seq}
    sequence_file = DATA_DIR / "sequence.fasta"
    sequence_file.parent.mkdir(parents=True, exist_ok=True)
    with sequence_file.open("w") as handle:
        for strain, sequence in seq.items():
            handle.write(f">{strain}\n{sequence}\n")

    print(len(seq.keys()), len(auspice["strain"]))
    print(set(seq.keys()) - set(auspice["strain"]))
    print(set(auspice["strain"]) - set(seq.keys()))

    auspice.to_csv(DATA_DIR / "auspice_metadata.tsv", sep="\t", index=False)

    consensus_seq = {}
    for sample, path in samples[["sample", "nextstrain_genome"]].loc[samples["nextstrain_genome"].notnull()].values:
        if path in {"unknown", "-"}:
            continue
        fasta_path = Path(path.replace("/analysis_MN90894/", "/analysis_MN908947/"))
        with fasta_path.open("r") as handle:
            key = ""
            for raw_line in handle:
                line = raw_line.strip()
                if line.startswith(">"):
                    key = sample
                    consensus_seq[key] = ""
                    continue
                consensus_seq[key] += line

    wave1_consensus = WAVE1_PATH / "1wave_all_consensus_sequences.fasta"
    wave1_data: Dict[str, str] = {}
    with wave1_consensus.open("r") as handle:
        key = ""
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                key = line[1:]
                wave1_data[key] = ""
                continue
            wave1_data[key] += line

    seq = {**consensus_seq, **wave1_data}
    output = DATA_DIR / "AND1-100_1wave_consensus_sequences.fasta"
    with output.open("w") as handle:
        for strain, sequence in seq.items():
            handle.write(f">{strain}\n{sequence}\n")


def main() -> None:
    plan = read_batch_plan()
    samples = build_samples(plan)
    lineage = build_lineages(plan)
    nextclade = build_nextclade(plan)
    metadata = build_metadata()

    samples = samples.loc[samples["sample"] != "AND37139"]
    samples = merge_metadata(samples, lineage, nextclade, metadata)
    samples = assign_geography(samples)
    samples = assign_hospitals(samples)
    samples = update_dates(samples)

    complete_path = SSP_PATH / "analysis_MN908947/other/complete_samples_metadata.tsv"
    complete_path.parent.mkdir(parents=True, exist_ok=True)
    samples.to_csv(complete_path, sep="\t")

    auspice = build_auspice(samples, nextclade, lineage)
    write_sequences(auspice, samples)


if __name__ == "__main__":
    main()
