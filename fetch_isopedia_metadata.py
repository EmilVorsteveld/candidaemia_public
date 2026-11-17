"""Fetch extra metadata for isopedia samples from SRA (EXPERIMENT_PACKAGE XML) and ENCODE.

Usage (from the project root):
  python3 -I scripts/fetch_isopedia_metadata.py \
    data/kinnex/sqanti/filter_out/sqanti_filter_out.filtered.tsv \
    data/ext_data/isopedia/isopedia_sra_encode_metadata.tsv

Reads the ##[SAMPLE] header of the isopedia output (first column = experiment accession
SRX/ERX/DRX/ENCSR, column 11 = source). Joined onto isopedia_samples in 05_transcriptome.Rmd.

single_cell = "yes" if SRA library_source is TRANSCRIPTOMIC SINGLE CELL, or the sample's own
title/protocol mentions single-cell/single-nuclei/10x Chromium (the "NEBNext Single Cell/Low
Input" bulk cDNA kit name is ignored).
"""
import csv, json, re, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA = {"User-Agent": "isopedia-metadata-fetch/0.1"}

def get(url, data=None, headers=None, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={**UA, **(headers or {})})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))

def text(el, path):
    x = el.find(path)
    return (x.text or "").strip() if x is not None and x.text else ""

SC_PAT = re.compile(r"single[- ]?cell|single[- ]?nucle|scRNA|snRNA|10x Genomics|chromium|scISO|cell barcode|_sc_|droplet", re.I)
LOW_INPUT = re.compile(r"single[- ]cell\s*/\s*low[- ]input", re.I)

def fetch_sra(accs, batch=150):
    out = {}
    for i in range(0, len(accs), batch):
        chunk = accs[i:i + batch]
        term = " OR ".join(f"{a}[ACCN]" for a in chunk)
        es = json.loads(get(EUTILS + "esearch.fcgi", data=urllib.parse.urlencode(
            {"db": "sra", "term": term, "usehistory": "y", "retmode": "json", "retmax": 10000}).encode()))["esearchresult"]
        xml = get(EUTILS + "efetch.fcgi", data=urllib.parse.urlencode(
            {"db": "sra", "WebEnv": es["webenv"], "query_key": es["querykey"], "retmode": "xml", "retmax": 10000}).encode())
        root = ET.fromstring(xml)
        for pkg in root.findall("EXPERIMENT_PACKAGE"):
            exp = pkg.find("EXPERIMENT")
            acc = exp.get("accession")
            ld = exp.find("DESIGN/LIBRARY_DESCRIPTOR")
            inst = exp.find(".//INSTRUMENT_MODEL")
            samp = pkg.find("SAMPLE")
            attrs = []
            if samp is not None:
                for a in samp.findall("SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
                    attrs.append((text(a, "TAG"), text(a, "VALUE")))
            runs = pkg.findall("RUN_SET/RUN")
            spots = sum(int(r.get("total_spots") or 0) for r in runs)
            bases = sum(int(r.get("total_bases") or 0) for r in runs)
            rec = {
                "experiment_title": text(exp, "TITLE"),
                "organism": text(samp, "SAMPLE_NAME/SCIENTIFIC_NAME") if samp is not None else "",
                "study_title": text(pkg, "STUDY/DESCRIPTOR/STUDY_TITLE"),
                "sample_title": text(samp, "TITLE") if samp is not None else "",
                "library_name": text(ld, "LIBRARY_NAME"),
                "library_strategy": text(ld, "LIBRARY_STRATEGY"),
                "library_source": text(ld, "LIBRARY_SOURCE"),
                "library_selection": text(ld, "LIBRARY_SELECTION"),
                "library_layout": (ld.find("LIBRARY_LAYOUT")[0].tag if ld is not None and ld.find("LIBRARY_LAYOUT") is not None and len(ld.find("LIBRARY_LAYOUT")) else ""),
                "library_protocol": text(ld, "LIBRARY_CONSTRUCTION_PROTOCOL"),
                "design_description": text(exp, "DESIGN/DESIGN_DESCRIPTION"),
                "instrument": (inst.text or "").strip() if inst is not None else "",
                "submitter": text(pkg, "Organization/Name"),
                "n_runs": len(runs),
                "run_accessions": ";".join(r.get("accession") for r in runs),
                "sra_spots": spots,
                "sra_bases": bases,
                "sample_attributes": "; ".join(f"{t}={v}" for t, v in attrs),
            }
            out[acc] = rec
        print(f"SRA {min(i + batch, len(accs))}/{len(accs)}", file=sys.stderr)
        time.sleep(0.5)
    return out

def fetch_encode(accs, batch=50):
    out = {}
    for i in range(0, len(accs), batch):
        q = [("type", "Experiment"), ("format", "json"), ("limit", "all"), ("frame", "embedded")]
        q += [("accession", a) for a in accs[i:i + batch]]
        d = json.loads(get("https://www.encodeproject.org/search/?" + urllib.parse.urlencode(q),
                           headers={"Accept": "application/json"}))
        for e in d.get("@graph", []):
            libs = [r.get("library") or {} for r in e.get("replicates", [])]
            bios = [l.get("biosample") or {} for l in libs]
            donors = [b.get("donor") or {} for b in bios]
            def uniq(xs):
                return ";".join(sorted({str(x) for x in xs if x not in (None, "", [])}))
            out[e["accession"]] = {
                "experiment_title": e.get("assay_title", ""),
                "organism": uniq(((d.get("organism") or {}).get("scientific_name") if isinstance(d.get("organism"), dict) else "Homo sapiens") for d in donors),
                "study_title": (e.get("lab") or {}).get("title", ""),
                "sample_title": e.get("biosample_summary", ""),
                "library_strategy": e.get("assay_term_name", ""),
                "library_source": uniq(l.get("nucleic_acid_term_name") for l in libs),
                "library_protocol": uniq(l.get("construction_method") if isinstance(l.get("construction_method"), str) else ",".join(l.get("construction_method") or []) for l in libs),
                "design_description": e.get("description", ""),
                "instrument": "",
                "submitter": (e.get("lab") or {}).get("title", ""),
                "n_runs": len(e.get("replicates", [])),
                "sample_attributes": "; ".join(filter(None, [
                    "biosample_type=" + (e.get("biosample_ontology") or {}).get("classification", ""),
                    "biosample=" + (e.get("biosample_ontology") or {}).get("term_name", ""),
                    "life_stage_age=" + str(e.get("life_stage_age", "")),
                    "sex=" + uniq(d.get("sex") for d in donors),
                    "age=" + uniq(d.get("age") for d in donors),
                    "treatments=" + uniq(t.get("treatment_term_name") for b in bios for t in (b.get("treatments") or [])),
                ])),
            }
        print(f"ENCODE {min(i + batch, len(accs))}/{len(accs)}", file=sys.stderr)
    return out

def main(inp, outp):
    rows = []
    with open(inp) as f:
        for line in f:
            if not line.startswith("##"):
                break
            if line.startswith("##[SAMPLE]"):
                rows.append(line[len("##[SAMPLE]"):].rstrip("\n").split("\t"))
    rows = rows[1:]  # drop the column-name line
    sra = [r[0] for r in rows if r[10] == "SRA"]
    enc = [r[0] for r in rows if r[10] == "ENCODE"]
    meta = {**fetch_sra(sra), **fetch_encode(enc)}
    cols = ["experiment_title", "organism", "study_title", "sample_title", "library_name", "library_strategy", "library_source",
            "library_selection", "library_layout", "library_protocol", "design_description", "instrument", "submitter",
            "n_runs", "run_accessions", "sra_spots", "sra_bases", "sample_attributes"]
    with open(outp, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "found"] + cols + ["single_cell"])
        for r in rows:
            m = meta.get(r[0], {})
            blob = LOW_INPUT.sub("", " ".join(str(m.get(c, "")) for c in ["experiment_title", "sample_title", "library_name", "library_protocol", "design_description"]))
            sc = "yes" if (m.get("library_source", "") == "TRANSCRIPTOMIC SINGLE CELL" or SC_PAT.search(blob)) else "no"
            w.writerow([r[0], "yes" if m else "no"] + [re.sub(r"\s+", " ", str(m.get(c, ""))) for c in cols] + [sc])
    print(f"found {len(meta)}/{len(rows)}", file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
