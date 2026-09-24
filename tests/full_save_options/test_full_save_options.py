"""
Checks the full save options set as attributes of <save><full_data><enable>:

    <enable microenvironment="false" cell_data="true" neighbor_data="false"
            attachments="true" string_attachments="true"
            intracellular_data="true">true</enable>

For each case, a config is generated from the template project, the driver
(built with `make`) writes two full saves, and we check that exactly the
enabled parts are written, both as files and as entries in the output XML.

intracellular_data is only written in PhysiBoSS builds, so here we only check
that turning it off does not affect the other parts.

If pcdl is installed, each output is also loaded with pcdl.TimeSeries, with
its load options matching what was saved:
    microenv  <- microenvironment
    graph     <- neighbor_data and attachments and string_attachments
    physiboss <- intracellular_data

Usage: python3 test_full_save_options.py [path/to/full_save_options]
"""

import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

try:
    import pcdl
except ImportError:
    pcdl = None

HERE = os.path.dirname(os.path.abspath(__file__))
BASE_CONFIG = os.path.join(HERE, "..", "..", "sample_projects", "template", "config", "PhysiCell_settings.xml")
DRIVER = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "full_save_options")

CUSTOM = "cellular_information/cell_populations/cell_population/custom"

# option name -> (file suffix, path of the <filename> node in the output XML)
PARTS = {
    "microenvironment": ("_microenvironment0.mat", "microenvironment/domain/data/filename"),
    "cell_data": ("_cells.mat", CUSTOM + "/simplified_data/filename"),
    "neighbor_data": ("_cell_neighbor_graph.txt", CUSTOM + "/neighbor_graph/filename"),
    "attachments": ("_attached_cells_graph.txt", CUSTOM + "/attached_cells_graph/filename"),
    "string_attachments": ("_spring_attached_cells_graph.txt", CUSTOM + "/spring_attached_cells_graph/filename"),
}
ALL_OPTIONS = list(PARTS) + ["intracellular_data"]
GRAPH_PARTS = {"neighbor_data", "attachments", "string_attachments"}

N_SAVES = 2


def make_config(folder, attributes):
    tree = ET.parse(BASE_CONFIG)
    root = tree.getroot()
    root.find("save/folder").text = folder
    enable = root.find("save/full_data/enable")
    enable.text = "true"
    for key, value in attributes.items():
        enable.set(key, value)
    # pcdl reads the cell types from the settings file in the output folder
    path = os.path.join(folder, "PhysiCell_settings.xml")
    tree.write(path)
    return path


def run_case(name, attributes, expected):
    """expected: set of option names from PARTS that should be saved."""
    errors = []
    with tempfile.TemporaryDirectory(prefix="full_save_") as folder:
        config = make_config(folder, attributes)
        result = subprocess.run([DRIVER, config], cwd=folder, capture_output=True, text=True)
        if result.returncode != 0:
            return ["driver exited with code %d:\n%s%s" % (result.returncode, result.stdout[-2000:], result.stderr[-2000:])]

        for n in range(N_SAVES):
            base = "output%08d" % n
            xml_path = os.path.join(folder, base + ".xml")
            if not os.path.exists(xml_path):
                errors.append("%s.xml was not written" % base)
                continue
            try:
                xml_root = ET.parse(xml_path).getroot()
            except ET.ParseError as e:
                errors.append("%s.xml does not parse: %s" % (base, e))
                continue

            for option, (suffix, xml_node) in PARTS.items():
                filename = base + suffix
                file_exists = os.path.exists(os.path.join(folder, filename))
                node = xml_root.find(xml_node)
                if option in expected:
                    if not file_exists:
                        errors.append("%s: missing %s" % (option, filename))
                    if node is None:
                        errors.append("%s: missing <%s> in %s.xml" % (option, xml_node, base))
                    elif (node.text or "").strip() != filename:
                        errors.append("%s: <%s> is '%s' in %s.xml, expected '%s'"
                                      % (option, xml_node, node.text, base, filename))
                else:
                    if file_exists:
                        errors.append("%s: %s was written but the option is off" % (option, filename))
                    if node is not None:
                        errors.append("%s: <%s> is in %s.xml but the option is off" % (option, xml_node, base))

            # the cell legend is part of the XML regardless of cell_data
            if xml_root.find(CUSTOM + "/simplified_data/labels") is None:
                errors.append("cell labels missing from %s.xml" % base)

            # the driver attaches cells 0-1, and spring-attaches cells 2-3
            for option, pair in (("attachments", {"0", "1"}), ("string_attachments", {"2", "3"})):
                if option not in expected:
                    continue
                with open(os.path.join(folder, base + PARTS[option][0])) as f:
                    ids = set()
                    for line in f:
                        line = line.strip()
                        if line and ":" in line:
                            cell, others = line.split(":", 1)
                            if others.strip():
                                ids.add(cell.strip())
                if ids != pair:
                    errors.append("%s: expected cells %s to have attachments in %s, got %s"
                                  % (option, sorted(pair), base, sorted(ids)))

            # the mesh is always saved, pcdl needs it even without the microenvironment
            if xml_root.find("microenvironment/domain/mesh/voxels/filename") is None:
                errors.append("mesh missing from %s.xml" % base)

        if pcdl is not None and "cell_data" in expected:
            errors += check_pcdl(folder, attributes, expected)
    return errors


def check_pcdl(folder, attributes, expected):
    """Load the output with pcdl, skipping the parts that were not saved."""
    microenv = "microenvironment" in expected
    graph = GRAPH_PARTS <= expected
    physiboss = attributes.get("intracellular_data", "true") != "false"
    try:
        mcdsts = pcdl.TimeSeries(folder, microenv=microenv, graph=graph, physiboss=physiboss, verbose=False)
        mcds_list = mcdsts.get_mcds_list()
    except Exception as e:
        return ["pcdl (microenv=%s, graph=%s, physiboss=%s) failed to load: %s: %s"
                % (microenv, graph, physiboss, type(e).__name__, e)]

    errors = []
    if len(mcds_list) != N_SAVES:
        errors.append("pcdl: loaded %d time steps, expected %d" % (len(mcds_list), N_SAVES))
    for mcds in mcds_list:
        n_cells = mcds.get_cell_df().shape[0]
        if n_cells != 4:
            errors.append("pcdl: %d cells at t=%s, expected 4" % (n_cells, mcds.get_time()))
        if microenv and "substrate" not in mcds.get_conc_df().columns:
            errors.append("pcdl: substrate concentrations missing at t=%s" % mcds.get_time())
        if graph:
            attached = {k: set(v) for k, v in mcds.get_attached_graph_dict().items() if len(v)}
            spring = {k: set(v) for k, v in mcds.get_spring_graph_dict().items() if len(v)}
            if attached != {0: {1}, 1: {0}}:
                errors.append("pcdl: attached graph is %s at t=%s" % (attached, mcds.get_time()))
            if spring != {2: {3}, 3: {2}}:
                errors.append("pcdl: spring graph is %s at t=%s" % (spring, mcds.get_time()))
    return errors


def main():
    if not os.path.exists(DRIVER):
        print("Driver not found at %s, run `make` first." % DRIVER)
        return 1
    if pcdl is None:
        print("pcdl is not installed, skipping the pcdl loading checks.\n")
    else:
        print("Also loading the outputs with pcdl %s.\n" % pcdl.__version__)

    cases = [("no attributes (defaults)", {}, set(PARTS)),
             ("all explicitly true", {k: "true" for k in ALL_OPTIONS}, set(PARTS)),
             ("all false", {k: "false" for k in ALL_OPTIONS}, set())]
    for option in ALL_OPTIONS:
        cases.append(("only %s off" % option, {option: "false"}, set(PARTS) - {option}))
        if option in PARTS:
            others_off = {k: "false" for k in ALL_OPTIONS if k != option}
            cases.append(("only %s on" % option, others_off, {option}))

    failed = 0
    for name, attributes, expected in cases:
        errors = run_case(name, attributes, expected)
        print("[%s] %s" % ("FAIL" if errors else " OK ", name))
        for e in errors:
            print("       " + e)
        failed += bool(errors)

    print("\n%d/%d cases passed" % (len(cases) - failed, len(cases)))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
