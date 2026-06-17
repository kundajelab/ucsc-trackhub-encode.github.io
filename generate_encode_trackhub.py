#!/usr/bin/env python3
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import urlparse


META_COLUMNS = {
    "experiment",
    "type",
    "assay",
    "target",
    "biosample",
    "tissue",
    "organ",
    "system",
    "model_annotation",
    "qc_flag",
}


def slug(value, fallback="unknown"):
    value = (value or "").strip()
    if value == "+":
        value = "plus"
    elif value == "-":
        value = "minus"
    value = re.sub(r"[^A-Za-z0-9]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or fallback


def track_name(*parts):
    name = "_".join(slug(p) for p in parts if p)
    if not re.match(r"^[A-Za-z]", name):
        name = "t_" + name
    return name[:240]


def label(value, max_len=17):
    value = re.sub(r"\s+", " ", (value or "").strip())
    return value[:max_len]


def pretty(value):
    value = (value or "").strip()
    return value if value else "unknown"


def subgroup_label(value):
    value = pretty(value).replace(" ", "_")
    value = re.sub(r"[^A-Za-z0-9_.+-]+", "_", value)
    return value or "unknown"


def clean_output_column(column):
    base = re.sub(r"\s+(bigWig|bigBed)\s*$", "", column, flags=re.I).strip()
    output = "NA"

    strand = re.search(r"\((plus|minus) strand\)", base, flags=re.I)
    if strand:
        output = strand.group(1).capitalize()
        base = re.sub(r"\s*\((plus|minus) strand\)", "", base, flags=re.I).strip()

    count_profile = re.match(r"^(counts|profile)\s+(.+)$", base, flags=re.I)
    if count_profile:
        output = count_profile.group(1).capitalize()
        base = count_profile.group(2).strip()

    return base, output


def file_type(url):
    path = urlparse(url).path.lower()
    if path.endswith(".bigwig") or path.endswith(".bw"):
        return "bigWig"
    if path.endswith(".bigbed") or path.endswith(".bb"):
        return "bigBed"
    return None


def accession_from_url(url):
    m = re.search(r"/files/([^/]+)/", url)
    if m:
        return m.group(1)
    stem = Path(urlparse(url).path).stem
    return stem or "file"


def subgroup_defs(records):
    fields = [
        ("model", "Model"),
        ("track", "Track"),
        ("assay", "Assay"),
        ("tissue", "Tissue"),
        ("organ", "Organ"),
        ("system", "System"),
        ("biosample", "Biosample"),
        ("output", "Output"),
        ("target", "Target"),
    ]
    lines = []
    for i, (field, title) in enumerate(fields, 1):
        vals = sorted({r[field] for r in records if r[field]})
        if not vals:
            vals = ["unknown"]
        pairs = " ".join(f"{slug(v)}={subgroup_label(v)}" for v in vals)
        lines.append(f"subGroup{i} {field} {title} {pairs}")
    return lines


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: generate_encode_trackhub.py metadata.tsv hub_repo_dir")

    metadata_path = Path(sys.argv[1])
    repo_dir = Path(sys.argv[2])
    hg_dir = repo_dir / "hg38"
    hg_dir.mkdir(parents=True, exist_ok=True)
    for old in hg_dir.glob("trackDb.*.txt"):
        old.unlink()

    with metadata_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    url_columns = [c for c in fieldnames if c not in META_COLUMNS]
    groups = defaultdict(list)

    for row in rows:
        model = pretty(row.get("type"))
        assay = pretty(row.get("assay"))
        target = pretty(row.get("target"))
        biosample = pretty(row.get("biosample"))
        tissue = pretty(row.get("tissue"))
        organ = pretty(row.get("organ"))
        system = pretty(row.get("system"))
        experiment = pretty(row.get("experiment"))
        model_annotation = pretty(row.get("model_annotation"))

        for output_col in url_columns:
            url = (row.get(output_col) or "").strip()
            if not url:
                continue
            ttype = file_type(url)
            if ttype not in {"bigWig", "bigBed"}:
                continue
            output_family, output = clean_output_column(output_col)
            accession = accession_from_url(url)
            record = {
                "track_id": track_name(model, output_family, output, accession),
                "type": ttype,
                "url": url,
                "accession": accession,
                "experiment": experiment,
                "model_annotation": model_annotation,
                "model": model,
                "assay": assay,
                "target": target,
                "biosample": biosample,
                "tissue": tissue,
                "organ": organ,
                "system": system,
                "output_family": output_family,
                "track": output_family,
                "output": output,
            }
            groups[(model, output_family, ttype)].append(record)

    includes = []
    for (model, output, ttype), records in sorted(groups.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
        comp_id = track_name(model, output, ttype)
        filename = f"trackDb.{slug(model).lower()}.{slug(output).lower()}.{ttype.lower()}.txt"
        includes.append(filename)
        comp_short = label(f"{model} {output}", 17)
        type_line = "type bigWig" if ttype == "bigWig" else "type bigBed 6"
        visibility = "full" if ttype == "bigWig" else "dense"
        color = "0,90,180" if ttype == "bigWig" else "180,60,0"

        out = []
        out.extend([
            f"track {comp_id}",
            "compositeTrack on",
            "parent EncodeDLTrack",
            f"shortLabel {comp_short}",
            f"longLabel {model} {output} {ttype} tracks",
            type_line,
            "visibility hide",
            "allButtonPair on",
            "dragAndDrop subTracks",
        ])
        out.extend(subgroup_defs(records))
        out.extend([
            "dimensions dimX=target dimY=biosample dimA=model dimB=assay",
            "filterComposite dimA dimB",
            "sortOrder model=+ assay=+ tissue=+ organ=+ system=+ biosample=+ output=+ target=+",
            "",
        ])

        seen = set()
        for r in sorted(records, key=lambda rec: (rec["target"], rec["biosample"], rec["assay"], rec["output"], rec["accession"])):
            tid = r["track_id"]
            if tid in seen:
                tid = track_name(tid, r["experiment"])
            seen.add(tid)
            short = label(f"{r['target']} {r['biosample']} {r['output']}", 17)
            long = f"{r['model']} {r['output_family']} {r['output']} {r['target']} {r['biosample']} {r['assay']} {r['accession']}"
            out.extend([
                f"track {tid}",
                f"parent {comp_id} off",
                type_line,
                f"bigDataUrl {r['url']}",
                f"shortLabel {short}",
                f"longLabel {long}",
                (
                    "subGroups "
                    f"model={slug(r['model'])} "
                    f"track={slug(r['track'])} "
                    f"assay={slug(r['assay'])} "
                    f"tissue={slug(r['tissue'])} "
                    f"organ={slug(r['organ'])} "
                    f"system={slug(r['system'])} "
                    f"biosample={slug(r['biosample'])} "
                    f"output={slug(r['output'])} "
                    f"target={slug(r['target'])}"
                ),
                f"visibility {visibility}",
                f"color {color}",
            ])
            if ttype == "bigWig":
                out.append("autoScale on")
            out.append("")

        (hg_dir / filename).write_text("\n".join(out), encoding="utf-8")

    (hg_dir / "trackDb.txt").write_text(
        "\n".join([
            "track EncodeDLTrack",
            "superTrack on show",
            "shortLabel ENCODE DL",
            "longLabel ENCODE deep learning sequence models: BPNet, ChromBPNet, ProCapNet, ReporterNet",
            "",
            *[f"include {name}" for name in includes],
            "",
        ]),
        encoding="utf-8",
    )

    (repo_dir / "hub.txt").write_text(
        "\n".join([
            "hub EncodeDLHub",
            "shortLabel ENCODE DL",
            "longLabel A collection of model predictions, contribution scores, and motif instances from ENCODE deep learning sequence models: BPNet, ChromBPNet, ProCapNet, ReporterNet",
            "genomesFile genomes.txt",
            "email akundaje@stanford.edu",
            "",
        ]),
        encoding="utf-8",
    )
    (repo_dir / "genomes.txt").write_text("genome hg38\ntrackDb hg38/trackDb.txt\n", encoding="utf-8")
    (repo_dir / ".nojekyll").write_text("", encoding="utf-8")

    print(f"wrote {len(groups)} composites")
    print(f"wrote {sum(len(v) for v in groups.values())} bigWig/bigBed subtracks")


if __name__ == "__main__":
    main()
