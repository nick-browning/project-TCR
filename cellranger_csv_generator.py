#!/usr/bin/env python3
"""Shiny app for building a cellranger multi config CSV.

Paste one SRA accession or a whole list (one per line, or separated by spaces
or commas), pick their data type from the dropdown, and the app
writes the [libraries] section of cellrangermulti_config.csv pointing at the
FASTQs that sra_tools_final.sh renamed.

Run it with:

    python3 cellranger_csv_generator.py

then open the URL it prints. An unused port is chosen automatically, because
this login node is shared and the usual ports are generally already serving
other people's work; pass --port to pick one yourself. On VS Code remote the
port is forwarded automatically, so the printed link works straight from the
terminal.

Sections whose libraries are absent are left out, because cellranger multi
rejects a [vdj] or [feature] reference with no matching library.

The file is in four parts:

    1. Settings      -- the paths and data types this project uses.
    2. Python functions  -- find the FASTQs, build the CSV text, look for mistakes.
    3. The interface -- what the page looks like, and what reacts to what.
    4. Startup       -- pick a port and serve the app.
"""

import argparse
import re
import socket
from pathlib import Path

import shiny
from shiny import App, reactive, render, ui


# --------------------------------------------------------------------------
# 1. Settings
# Defaults for this project, paths can be edited in the shiny interface. 
# Editing a box in the shiny interface changes the path for the current session, but below are the default values.
# --------------------------------------------------------------------------

FASTQ_DIR = "/data/scratch/projects/punim1612/nick/sra_tools_sample_run"
OUTPUT_CSV = "/data/projects/punim1612/CSI/Nick-project/Slurm_scripts/cellrangermulti_config.csv"

GEX_REFERENCE = "/data/projects/punim1612/references/human/refdata-gex-GRCh38-2020-A"
VDJ_REFERENCE = "/data/projects/punim1612/references/human/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"
FEATURE_REFERENCE = "/data/projects/punim1612/CSI/Nick-project/cellranger_multi_feature_ref.csv"

# What you pick in the dropdown, the feature_types value cellranger wants.
# CITE-seq and generic feature barcoding are both antibody capture as far as
# cellranger is concerned; BEAM is antigen capture.

# TODO - need to be able to select multiple data types for a single accession, e.g. GEX + CITE-seq. 
# Currently the dropdown only allows one selection per accession.
DATA_TYPES = {
    "GEX": "Gene Expression",
    "TCR": "VDJ-T",
    "BCR": "VDJ-B",
    "CITE": "Antibody Capture",
    "FB": "Antibody Capture",
    "BEAM": "Antigen Capture",
}

# The dropdown shows both labels in DATA_TYPES for clarity. 
DROPDOWN_LABELS = {key: f"{key} - {value}" for key, value in DATA_TYPES.items()}

# Which feature_types values require which config section. 
# VDJ section is needed for VDJ-T and VDJ-B, feature section is needed for Antibody Capture and Antigen Capture.
VDJ_TYPES = {"VDJ-T", "VDJ-B"}
FEATURE_TYPES = {"Antibody Capture", "Antigen Capture"}

# An SRA / ENA / DDBJ run accession: SRR, ERR or DRR followed by digits.
# Anything else typed into the accession box is rejected rather than added,
# e.g. two accessions run together as SRR35210236SRR35210235.
ACCESSION_PATTERN = re.compile(r"^[SED]RR\d+$")


# --------------------------------------------------------------------------
# 2. Python functions
# Look for an accession's FASTQs on disk, turn a list of libraries
# into CSV text, and list anything that looks wrong. A "library" throughout is
# just a dict: {"accession": "SRR9729407", "type": "TCR"}.
# --------------------------------------------------------------------------

def find_fastqs(fastq_dir, accession):
    """Return the reads present for an accession, e.g. {"R1", "R2"}."""
    directory = Path(fastq_dir)
    if not directory.is_dir():
        return set()  #Checks if fastq directory exists, if not returns empty set (skips rest of function)

    reads = set()  #otherwise if non-empty directory, looks for fastq files 
                   #matching the pattern <accession>_S*_L*_*_001.fastq.gz and adds the read (R1 or R2) to the set of reads found.
    for path in directory.glob(f"{accession}_S*_L*_*_001.fastq.gz"):
        # Loops over all files in the directory that match the pattern for the given accession.
        # Splits the filename by underscores and takes the second-last piece, which should be "R1" or "R2" (-1 is 001.fastq.gz, -2 is R1 or R2).
        reads.add(path.name.split("_")[-2])
    return reads


def build_config(libraries, fastq_dir, gex_ref, vdj_ref, feature_ref):
    """Render the cellranger multi config CSV for the given libraries.

    Each section is emitted only when a library needs it -- cellranger multi
    errors out on a [vdj] or [feature] reference with nothing to use it.
    """
    features = {DATA_TYPES[lib["type"]] for lib in libraries}
    lines = []

    if "Gene Expression" in features:
        lines += ["[gene-expression]", f"reference,{gex_ref}", "create-bam,true", ""]

    if features & VDJ_TYPES:
        lines += ["[vdj]", f"reference,{vdj_ref}", ""]

    if features & FEATURE_TYPES:
        lines += ["[feature]", f"reference,{feature_ref}", ""]

    lines += ["[libraries]", "fastq_id,fastqs,feature_types"]
    for lib in libraries:
        lines.append(f"{lib['accession']},{fastq_dir},{DATA_TYPES[lib['type']]}")

    return "\n".join(lines) + "\n"


def find_problems(libraries, fastq_dir, gex_ref, vdj_ref, feature_ref):
    """Return (errors, warnings). Errors block writing; warnings just show."""
    errors = []
    warnings = []

    if not libraries:
        return ["Add at least one library."], []

    accessions = [lib["accession"] for lib in libraries]
    features = {DATA_TYPES[lib["type"]] for lib in libraries}

    # Does every accession have both of its reads on disk?
    if not Path(fastq_dir).is_dir():
        errors.append(f"FASTQ directory does not exist: {fastq_dir}")
    else:
        for accession in dict.fromkeys(accessions):   # each one once
            reads = find_fastqs(fastq_dir, accession)
            if not reads:
                errors.append(
                    f"No FASTQs found for {accession}. Expected "
                    f"{accession}_S1_L001_R1_001.fastq.gz in the FASTQ directory."
                )
            elif not {"R1", "R2"} <= reads:
                missing = " and ".join(sorted({"R1", "R2"} - reads))
                errors.append(
                    f"{accession} is missing {missing}. Check the "
                    "sra_tools_final.sh .out/.err log: it leaves files unrenamed "
                    "if it doesn't find exactly two reads of 20 bp or longer."
                )

    # cellranger needs unique fastq_ids.
    for accession in sorted({a for a in accessions if accessions.count(a) > 1}):
        errors.append(f"{accession} is listed more than once.")

    # Combinations cellranger multi will not make sense of.
    if features & FEATURE_TYPES and "Gene Expression" not in features:
        warnings.append(
            "CITE / FB / BEAM libraries need a Gene Expression library in the "
            "same cellranger multi run."
        )
    if "Antigen Capture" in features and not features & VDJ_TYPES:
        warnings.append(
            "BEAM needs a TCR library in the same run to pair antigen barcodes "
            "with clonotypes."
        )

    # Only the references the chosen libraries actually point at.
    needed = []
    if "Gene Expression" in features:
        needed.append(("transcriptome reference", gex_ref))
    if features & VDJ_TYPES:
        needed.append(("VDJ reference", vdj_ref))
    if features & FEATURE_TYPES:
        needed.append(("feature reference", feature_ref))

    for name, path in needed:
        if not path:
            errors.append(f"No {name} given.")
        elif not Path(path).exists():
            hint = ""
            if name == "feature reference":
                hint = " Run cellranger_multi_feature_antibody_generator.py to create it."
            warnings.append(f"{name} not found: {path}.{hint}")

    return errors, warnings


def describe_reads(fastq_dir, accession):
    """Return (text, bootstrap colour class) for the table's FASTQs column."""
    reads = find_fastqs(fastq_dir, accession)
    found = ", ".join(sorted(reads))
    if {"R1", "R2"} <= reads:
        return found, "text-success"
    if reads:
        return f"{found} (incomplete)", "text-danger"
    return "none found", "text-danger"


def count_libraries(count):
    """1 -> "1 library", 2 -> "2 libraries"."""
    return "1 library" if count == 1 else f"{count} libraries"


# --------------------------------------------------------------------------
# 3. The interface
#
# A Shiny app is two halves that talk to each other through names:
#
#   app_ui  -- the page. Every ui.input_*/ui.output_* takes an id as its first
#              argument, and that id is the name the other half uses:
#              ui.input_text("accession", ...) is read as input.accession(),
#              and ui.output_ui("preview") is filled by the function below
#              called preview().
#   server  -- the behaviour, written as small functions that Shiny calls for
#              you whenever something they depend on changes.
#
# The class_= strings ("btn-primary", "alert alert-danger") are Bootstrap CSS
# class names. They only affect colour and spacing; nothing breaks if you
# change or drop them.
# --------------------------------------------------------------------------

def alert(text, kind):
    """A coloured message box. kind is success / danger / warning."""
    return ui.div(text, class_=f"alert alert-{kind} py-2")


def alert_list(heading, messages, kind):
    """A coloured box holding a bulleted list of messages."""
    return ui.div(
        ui.tags.strong(heading),
        ui.tags.ul(*[ui.tags.li(message) for message in messages]),
        class_=f"alert alert-{kind} py-2",
    )


app_ui = ui.page_sidebar(
    # Everything you click or type lives in the sidebar on the left.
    ui.sidebar(
        ui.h5("Add a library"),
        # A text area rather than input_text: a single-line box silently drops
        # the line breaks from a pasted list, gluing the accessions together.
        ui.input_text_area(
            "accession", "SRA accession(s)",
            placeholder="SRR35210236\nSRR35210235", rows=6,
        ),
        ui.input_select("data_type", "Data type", choices=DROPDOWN_LABELS, selected="GEX"),
        ui.input_action_button("add", "Add libraries", class_="btn-primary w-100"),
        ui.help_text(
            "Paste a list of accessions, one per line or separated by spaces or "
            "commas. They are all added with the selected data type."
        ),
        ui.hr(),
        ui.input_select("to_remove", "Remove libraries", choices=[], multiple=True),
        ui.input_action_button("remove", "Remove selected", class_="btn-outline-danger w-100"),
        ui.hr(),
        ui.h5("Paths"),
        ui.input_text("fastq_dir", "FASTQ directory", value=FASTQ_DIR),
        ui.input_text("gex_ref", "Transcriptome reference", value=GEX_REFERENCE),
        ui.input_text("vdj_ref", "VDJ reference", value=VDJ_REFERENCE),
        ui.input_text("feature_ref", "Feature reference", value=FEATURE_REFERENCE),
        ui.input_text("output", "Write config to", value=OUTPUT_CSV),
        ui.input_action_button("write", "Write config", class_="btn-success w-100"),
        width=380,
    ),
    # The main panel is read-only: three areas filled in by the server half.
    ui.h4("cellranger multi config generator"),
    ui.output_ui("status"),
    ui.h5("Libraries"),
    ui.output_ui("library_table"),
    ui.h5("cellrangermulti_config.csv preview"),
    ui.output_ui("preview"),
    title="cellranger multi config generator",
    fillable=False,
)


def server(input, output, session):
    # Four kinds of function appear below, and Shiny decides when to run each:
    #
    #   reactive.value   remembered state. Read it by calling it, libraries();
    #                    change it with libraries.set(...), and everything that
    #                    read it redraws itself.
    #   @reactive.effect does something (writes a file, shows a message). With
    #                    @reactive.event(input.add) it runs on that click only.
    #   @reactive.calc   works a value out from other values. Cached, so the
    #                    checks below run once per change, not once per reader.
    #   @render.ui       returns the HTML for the ui.output_ui of the same name.

    libraries = reactive.value([])   # list of {"accession", "type"}
    written = reactive.value("")     # status line from the last successful write

    # ---- adding and removing -------------------------------------------

    # ignore_init keeps these from firing on page load: an action button starts
    # at 0 rather than None, so without it opening the app would add a library
    # and overwrite the config straight away.
    @reactive.effect
    @reactive.event(input.add, ignore_init=True)
    def _add():
        # Newlines, tabs (a pasted spreadsheet column), spaces, commas and
        # semicolons all separate accessions.
        typed = [a.upper() for a in re.split(r"[\s,;]+", input.accession() or "") if a]
        if not typed:
            ui.notification_show("Type an SRA accession first.", type="warning")
            return

        invalid = list(dict.fromkeys(a for a in typed if not ACCESSION_PATTERN.match(a)))
        typed = [a for a in typed if ACCESSION_PATTERN.match(a)]

        current = list(libraries())
        already = {lib["accession"] for lib in current}
        # dict.fromkeys drops repeats within the pasted text, keeping the order.
        new = list(dict.fromkeys(a for a in typed if a not in already))
        skipped = list(dict.fromkeys(a for a in typed if a in already))

        current += [{"accession": a, "type": input.data_type()} for a in new]
        libraries.set(current)
        written.set("")
        ui.update_text_area("accession", value="")

        if invalid:
            ui.notification_show(
                f"Not an SRA run accession, not added: {', '.join(invalid)}",
                type="warning", duration=10,
            )
        if skipped:
            ui.notification_show(f"Already added: {', '.join(skipped)}", type="warning")
        if new:
            ui.notification_show(
                f"Added {count_libraries(len(new))} as {input.data_type()}.",
                type="message",
            )

    @reactive.effect
    @reactive.event(input.remove, ignore_init=True)
    def _remove():
        drop = set(input.to_remove() or [])
        if not drop:
            ui.notification_show("Select libraries to remove.", type="warning")
            return
        libraries.set([lib for lib in libraries() if lib["accession"] not in drop])
        written.set("")

    # Keep the "Remove libraries" list in step with what has been added.
    @reactive.effect
    def _sync_remove_choices():
        ui.update_select("to_remove", choices=[lib["accession"] for lib in libraries()])

    # ---- derived state -------------------------------------------------

    @reactive.calc
    def paths():
        """The path boxes, trimmed, ready to hand to the functions in part 2."""
        return {
            "fastq_dir": input.fastq_dir().strip(),
            "gex_ref": input.gex_ref().strip(),
            "vdj_ref": input.vdj_ref().strip(),
            "feature_ref": input.feature_ref().strip(),
        }

    @reactive.calc
    def problems():
        return find_problems(libraries(), **paths())

    @reactive.calc
    def config_text():
        if not libraries():
            return ""
        return build_config(libraries(), **paths())

    # ---- writing -------------------------------------------------------

    @reactive.effect
    @reactive.event(input.write, ignore_init=True)
    def _write():
        errors, _ = problems()
        if errors:
            ui.notification_show(
                "Fix the errors above before writing.", type="error", duration=8
            )
            return

        target = Path(input.output().strip())
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(config_text())
        except OSError as exc:
            written.set("")
            ui.notification_show(f"Could not write {target}: {exc}", type="error")
            return

        written.set(f"Wrote {target} with {count_libraries(len(libraries()))}.")
        ui.notification_show(f"Wrote {target}", type="message", duration=6)

    # ---- the three output areas -----------------------------------------

    @render.ui
    def status():
        """Any errors and warnings, plus the result of the last write."""
        errors, warnings = problems()
        boxes = []

        if written():
            boxes.append(alert(written(), "success"))
        if errors:
            boxes.append(alert_list("Errors", errors, "danger"))
        if warnings:
            boxes.append(alert_list("Warnings", warnings, "warning"))
        if libraries() and not errors and not warnings:
            boxes.append(alert("Ready to write.", "success"))

        return ui.TagList(*boxes)

    @render.ui
    def library_table():
        """One row per library, with whether its FASTQs were found."""
        if not libraries():
            return ui.p("No libraries yet.", class_="text-muted")

        rows = []
        for lib in libraries():
            found, colour = describe_reads(paths()["fastq_dir"], lib["accession"])
            rows.append(ui.tags.tr(
                ui.tags.td(ui.tags.code(lib["accession"])),
                ui.tags.td(lib["type"]),
                ui.tags.td(DATA_TYPES[lib["type"]]),
                ui.tags.td(found, class_=colour),
            ))

        return ui.tags.table(
            ui.tags.thead(ui.tags.tr(
                ui.tags.th("fastq_id"),
                ui.tags.th("Data type"),
                ui.tags.th("feature_types"),
                ui.tags.th("FASTQs"),
            )),
            ui.tags.tbody(*rows),
            class_="table table-sm table-striped",
        )

    @render.ui
    def preview():
        """The CSV exactly as the Write config button would save it."""
        text = config_text()
        if not text:
            return ui.p("Add a library to see the config.", class_="text-muted")
        return ui.tags.pre(
            text,
            class_="border rounded p-3 bg-light",
            style="white-space: pre; overflow-x: auto;",
        )


app = App(app_ui, server)


# --------------------------------------------------------------------------
# 4. Startup
# --------------------------------------------------------------------------

def free_port():
    """Ask the kernel for an unused port.

    The login node is shared, so the common ports are usually already taken by
    other people's servers -- binding a fixed port lands you on someone else's
    web root rather than this app.
    """
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        return probe.getsockname()[1]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--port", type=int, default=None,
        help="port to serve on (default: an unused one chosen automatically)",
    )
    args = parser.parse_args()

    port = args.port or free_port()
    # flush, or uvicorn's unbuffered logging prints over the top of the URL.
    print(
        "\n  cellranger multi config generator\n"
        f"  open http://127.0.0.1:{port}  (VS Code forwards this automatically)\n"
        "  press Ctrl-C to stop\n",
        flush=True,
    )

    shiny.run_app(app, host="127.0.0.1", port=port, launch_browser=False)
