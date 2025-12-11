#!/usr/bin/env python3
"""
dfba_configurator.py

Usage example:
    python dfba_configurator.py \
      --template config_template.xml \
      --output config_updated.xml \
      --config dfba_config.yaml \
      --sbml-folder ./models \
      --template-cell default_cell_def_name \
      --cell-prefix dfba_

Requirements:
    pip install lxml pandas pyyaml
    (libsbml is optional; not required for this version)
"""

import argparse
import os
import sys
import math
import pandas as pd
import glob
from lxml import etree
import copy
import yaml

# ----------------------
# Utilities
# ----------------------
def xml_pretty_write(tree, path):
    tree.write(path, pretty_print=True, xml_declaration=True, encoding="UTF-8")


def filename_stem(path):
    return os.path.splitext(os.path.basename(path))[0]


def find_template_cell(root, template_name=None):
    """
    Find a cell_definition by name, or return the first one if template_name is None or empty.
    
    Args:
        root: XML root element
        template_name: Optional name of the cell_definition. If None or empty string, returns first cell_definition.
    
    Returns:
        First matching cell_definition element
    
    Raises:
        ValueError: If template_name is provided but not found, or if no cell_definitions exist
    """
    if not template_name:
        # Return first cell_definition
        cell_def = root.find(".//cell_definition")
        if cell_def is None:
            raise ValueError("No cell_definition found in template.")
        return cell_def
    else:
        # Find by name
        cell_def = root.find(f".//cell_definition[@name='{template_name}']")
        if cell_def is None:
            raise ValueError(f"Template cell_definition with name '{template_name}' not found in template.")
        return cell_def


# ----------------------
# Microenvironment utilities
# ----------------------
def make_microenv_variable(name, diff=50000.0, decay=0.0, init=0.0):
    var = etree.Element("variable", name=name, units="dimensionless")
    var.set("ID", "-1")
    phys = etree.SubElement(var, "physical_parameter_set")
    etree.SubElement(phys, "diffusion_coefficient", units="micron^2/min").text = str(diff)
    etree.SubElement(phys, "decay_rate", units="1/min").text = str(decay)
    etree.SubElement(var, "initial_condition", units="mM").text = str(init)
    dbc = etree.SubElement(var, "Dirichlet_boundary_condition", units="mM", enabled="False")
    dbc.text = "0.0"
    do = etree.SubElement(var, "Dirichlet_options")
    for b in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"):
        etree.SubElement(do, "boundary_value", ID=b, enabled="False").text = "0.0"
    return var


# ----------------------
# Secretion block utilities (per cell_definition)
# ----------------------
def ensure_secretion_has_substrate(cell_def, substrate_name):
    """
    Inside cell_def -> phenotype -> secretion, ensure a <substrate name="..."> block exists.
    If not, create one with default zeros.
    """
    phenotype = cell_def.find("./phenotype")
    if phenotype is None:
        phenotype = etree.SubElement(cell_def, "phenotype")
    secretion = phenotype.find("secretion")
    if secretion is None:
        secretion = etree.SubElement(phenotype, "secretion")

    # search for substrate block with matching name attribute
    for sub in secretion.findall("substrate"):
        if sub.get("name") == substrate_name:
            # ensure fields exist
            _ensure_secretion_fields(sub)
            return

    # not found -> create default substrate block
    sub = etree.SubElement(secretion, "substrate", name=substrate_name)
    etree.SubElement(sub, "secretion_rate", units="1/min").text = "0.0"
    etree.SubElement(sub, "secretion_target", units="substrate density").text = "0.0"
    etree.SubElement(sub, "uptake_rate", units="1/min").text = "0.0"
    etree.SubElement(sub, "net_export_rate", units="total substrate/min").text = "0.0"


def _ensure_secretion_fields(sub_node):
    # Add any missing child tags with safe default values
    fields = {
        "secretion_rate": ("1/min", "0.0"),
        "secretion_target": ("substrate density", "0.0"),
        "uptake_rate": ("1/min", "0.0"),
        "net_export_rate": ("total substrate/min", "0.0")
    }
    for tag, (units, default) in fields.items():
        el = sub_node.find(tag)
        if el is None:
            etree.SubElement(sub_node, tag, units=units).text = default


# ----------------------
# Chemotactic sensitivities utilities (per cell_definition)
# ----------------------
def ensure_chemotactic_sensitivities(cell_def, substrate_name):
    """
    Ensure phenotype->motility->options->advanced_chemotaxis->chemotactic_sensitivities
    contains <chemotactic_sensitivity substrate="..." >0.0</chemotactic_sensitivity>
    """
    phenotype = cell_def.find("./phenotype")
    if phenotype is None:
        phenotype = etree.SubElement(cell_def, "phenotype")

    motility = phenotype.find("motility")
    if motility is None:
        motility = etree.SubElement(phenotype, "motility")

    options = motility.find("options")
    if options is None:
        options = etree.SubElement(motility, "options")

    adv = options.find("advanced_chemotaxis")
    if adv is None:
        adv = etree.SubElement(options, "advanced_chemotaxis")
        etree.SubElement(adv, "enabled").text = "false"
        etree.SubElement(adv, "normalize_each_gradient").text = "false"

    chems = adv.find("chemotactic_sensitivities")
    if chems is None:
        chems = etree.SubElement(adv, "chemotactic_sensitivities")

    # See if substrate already present
    for cs in chems.findall("chemotactic_sensitivity"):
        if cs.get("substrate") == substrate_name:
            # ensure text exists
            if (cs.text is None) or (cs.text.strip() == ""):
                cs.text = "0.0"
            return
    # Add missing
    new_cs = etree.SubElement(chems, "chemotactic_sensitivity", substrate=substrate_name)
    new_cs.text = "0.0"


# ----------------------
# Volume update utility
# ----------------------
def set_cell_volume(cell_def, total_volume=None, nuclear_volume=None):
    """
    Ensure phenotype->volume exists and update total and nuclear volumes if provided.
    """
    phenotype = cell_def.find("./phenotype")
    if phenotype is None:
        phenotype = etree.SubElement(cell_def, "phenotype")

    volume = phenotype.find("volume")
    if volume is None:
        volume = etree.SubElement(phenotype, "volume")

    if total_volume is not None:
        total_node = volume.find("total")
        if total_node is None:
            total_node = etree.SubElement(volume, "total", units="micron^3")
        total_node.text = str(total_volume)

    if nuclear_volume is not None:
        nuclear_node = volume.find("nuclear")
        if nuclear_node is None:
            nuclear_node = etree.SubElement(volume, "nuclear", units="micron^3")
        nuclear_node.text = str(nuclear_volume)


# ----------------------
# Intracellular dfba insertion
# ----------------------
def make_exchange_block_xml(exchange_config):
    """Create an exchange XML block from a YAML exchange configuration."""
    ex = etree.Element("exchange", substrate=str(exchange_config["substrate"]))
    etree.SubElement(ex, "fba_flux").text = str(exchange_config["fba_flux"])
    etree.SubElement(ex, "Km", units="mM").text = str(exchange_config["Km"])
    etree.SubElement(ex, "Vmax", units="fmol/pg DW cell/min").text = str(exchange_config["Vmax"])
    return ex


def add_intracellular_dfba(cell_def, sbml_path, model_config):
    """
    Adds (or replaces) phenotype->intracellular type="dfba" block in `cell_def`.
    
    Args:
        cell_def: XML element for cell_definition
        sbml_path: Path to SBML file
        model_config: Dictionary with model configuration from YAML:
            - exchanges: list of exchange dictionaries
            - growth_model: dict with cell_density, reference_volume, max_growth_rate, objective_reaction
            - settings: dict with intracellular_dt
            - death_model: dict with enabled flag
    """
    phenotype = cell_def.find("./phenotype")
    if phenotype is None:
        phenotype = etree.SubElement(cell_def, "phenotype")

    # Remove existing intracellular dfba node if any
    for old in phenotype.findall("intracellular"):
        if old.get("type") == "dfba":
            phenotype.remove(old)

    intracellular = etree.SubElement(phenotype, "intracellular", type="dfba")
    sbml_node = etree.SubElement(intracellular, "sbml_filename")
    sbml_node.text = sbml_path

    # Settings
    settings = etree.SubElement(intracellular, "settings")
    dt = model_config.get("settings", {}).get("intracellular_dt") if model_config.get("settings") else None
    if dt is None:
        dt = "0.01"
    etree.SubElement(settings, "intracellular_dt", units="min").text = str(dt)

    # Transport model (exchanges)
    transport = etree.SubElement(intracellular, "transport_model")
    exchanges = model_config.get("exchanges", [])
    for exchange in exchanges:
        transport.append(make_exchange_block_xml(exchange))

    # Growth model
    growth = etree.SubElement(intracellular, "growth_model")
    growth_config = model_config.get("growth_model", {})
    
    cell_density = growth_config.get("cell_density", "1.04")
    reference_volume = growth_config.get("reference_volume", "2494")
    nuclear_volume = growth_config.get("nuclear_volume", 0.0)
    max_growth_rate = growth_config.get("max_growth_rate", "0.86")
    objective_reaction = growth_config.get("objective_reaction", "R_Biomass")
    
    etree.SubElement(growth, "cell_density", units="g/ml").text = str(cell_density)
    etree.SubElement(growth, "reference_volume", units="pg").text = str(reference_volume)
    etree.SubElement(growth, "max_growth_rate", units="1/min").text = str(max_growth_rate)
    etree.SubElement(growth, "objective_reaction").text = str(objective_reaction)

    # Death model
    death_config = model_config.get("death_model", {})
    enabled = death_config.get("enabled", False) if death_config else False
    death_model = etree.SubElement(intracellular, "death_model", enabled=str(enabled).lower())
    death_type = death_config.get("death_type", "necrosis")
    death_trigger_flux = death_config.get("death_trigger_flux", "0.0")
    death_flux_threshold = death_config.get("death_flux_threshold", "0.0")
    death_rate_increase = death_config.get("death_rate_increase", "1.67e-5")

    etree.SubElement(death_model, "death_type", units="g/ml").text = str(death_type)
    etree.SubElement(death_model, "death_trigger_flux", units="1/min").text = str(death_trigger_flux)
    etree.SubElement(death_model, "death_flux_threshold").text = str(death_flux_threshold)
    etree.SubElement(death_model, "death_rate_increase").text = str(death_rate_increase)


# ----------------------
# Cell interactions matrix update
# ----------------------
def ensure_cell_interactions_for_all(cell_def, all_cell_names):
    """
    For a single cell_definition, ensure its cell_interactions section contains
    phagocytosis/attack/fusion/transformation entries for every name in all_cell_names.
    Missing entries are added with 0.0.
    """
    ci = cell_def.find("cell_interactions")
    if ci is None:
        ci = etree.SubElement(cell_def, "cell_interactions")
        # create some standard fields with defaults
        etree.SubElement(ci, "apoptotic_phagocytosis_rate", units="1/min").text = "0.0"
        etree.SubElement(ci, "necrotic_phagocytosis_rate", units="1/min").text = "0.0"
        etree.SubElement(ci, "other_dead_phagocytosis_rate", units="1/min").text = "0.0"

    # live_phagocytosis_rates -> phagocytosis_rate name="X"
    lpr = ci.find("live_phagocytosis_rates")
    if lpr is None:
        lpr = etree.SubElement(ci, "live_phagocytosis_rates")
    existing = {n.get("name") for n in lpr.findall("phagocytosis_rate")}
    for name in all_cell_names:
        if name not in existing:
            etree.SubElement(lpr, "phagocytosis_rate", name=name, units="1/min").text = "0.0"

    # attack_rates -> attack_rate name="X"
    ar = ci.find("attack_rates")
    if ar is None:
        ar = etree.SubElement(ci, "attack_rates")
    existing_attack = {n.get("name") for n in ar.findall("attack_rate")}
    for name in all_cell_names:
        if name not in existing_attack:
            etree.SubElement(ar, "attack_rate", name=name, units="1/min").text = "0.0"

    # fusion_rates -> fusion_rate name="X"
    fr = ci.find("fusion_rates")
    if fr is None:
        fr = etree.SubElement(ci, "fusion_rates")
    existing_fusion = {n.get("name") for n in fr.findall("fusion_rate")}
    for name in all_cell_names:
        if name not in existing_fusion:
            etree.SubElement(fr, "fusion_rate", name=name, units="1/min").text = "0.0"

    # transformation_rates -> transformation_rate name="X"
    ct = cell_def.find("cell_transformations")
    if ct is None:
        ct = etree.SubElement(cell_def, "cell_transformations")
    tr = ct.find("transformation_rates")
    if tr is None:
        tr = etree.SubElement(ct, "transformation_rates")
    existing_tr = {n.get("name") for n in tr.findall("transformation_rate")}
    for name in all_cell_names:
        if name not in existing_tr:
            etree.SubElement(tr, "transformation_rate", name=name, units="1/min").text = "0.0"


# ----------------------
# Validation
# ----------------------
def validate_yaml_config(config, config_yaml_path, sbml_folder):
    """
    Validate YAML configuration structure and content.
    
    Args:
        config: Parsed YAML configuration dictionary
        config_yaml_path: Path to YAML file (for error messages)
        sbml_folder: Folder containing SBML files (for path resolution)
    
    Raises:
        ValueError: If validation fails
    """
    errors = []
    warnings = []
    
    # Check required top-level sections
    if not isinstance(config, dict):
        raise ValueError(f"YAML config must be a dictionary, got {type(config)}")
    
    if "models" not in config:
        raise ValueError("YAML config must contain a 'models' section")
    
    if not isinstance(config["models"], dict):
        raise ValueError("'models' section must be a dictionary")
    
    if len(config["models"]) == 0:
        raise ValueError("'models' section must contain at least one model")
    
    # Build substrate name set
    substrate_names = set()
    if "substrates" in config:
        if not isinstance(config["substrates"], list):
            errors.append("'substrates' section must be a list")
        else:
            for idx, sub_config in enumerate(config["substrates"]):
                if not isinstance(sub_config, dict):
                    errors.append(f"Substrate entry {idx} must be a dictionary")
                    continue
                
                name = sub_config.get("name")
                if not name:
                    errors.append(f"Substrate entry {idx} missing required 'name' field")
                else:
                    if name in substrate_names:
                        warnings.append(f"Duplicate substrate name '{name}' in substrates section")
                    substrate_names.add(name)
                
                # Validate substrate properties
                if "diffusion_coefficient" in sub_config:
                    try:
                        diff = float(sub_config["diffusion_coefficient"])
                        if diff < 0:
                            errors.append(f"Substrate '{name}': diffusion_coefficient must be >= 0")
                    except (ValueError, TypeError):
                        errors.append(f"Substrate '{name}': diffusion_coefficient must be a number")
                
                if "decay_rate" in sub_config:
                    try:
                        decay = float(sub_config["decay_rate"])
                        if decay < 0:
                            errors.append(f"Substrate '{name}': decay_rate must be >= 0")
                    except (ValueError, TypeError):
                        errors.append(f"Substrate '{name}': decay_rate must be a number")
                
                if "initial_condition" in sub_config:
                    try:
                        float(sub_config["initial_condition"])
                    except (ValueError, TypeError):
                        errors.append(f"Substrate '{name}': initial_condition must be a number")
    
    # Validate models
    all_exchange_substrates = set()
    for model_name, model_config in config["models"].items():
        if not isinstance(model_config, dict):
            errors.append(f"Model '{model_name}': configuration must be a dictionary")
            continue
        
        # Check required fields
        if "sbml_path" not in model_config:
            errors.append(f"Model '{model_name}': missing required 'sbml_path' field")
        else:
            sbml_path_rel = model_config["sbml_path"]
            if not isinstance(sbml_path_rel, str):
                errors.append(f"Model '{model_name}': 'sbml_path' must be a string")
            else:
                # Resolve SBML path
                if os.path.isabs(sbml_path_rel):
                    sbml_path = sbml_path_rel
                else:
                    sbml_path = os.path.join(sbml_folder, sbml_path_rel)
                
                # Check if file exists
                if not os.path.exists(sbml_path):
                    errors.append(f"Model '{model_name}': SBML file not found: {sbml_path}")
                elif not os.path.isfile(sbml_path):
                    errors.append(f"Model '{model_name}': SBML path is not a file: {sbml_path}")
        
        # Validate exchanges
        if "exchanges" not in model_config:
            errors.append(f"Model '{model_name}': missing required 'exchanges' field")
        else:
            exchanges = model_config["exchanges"]
            if not isinstance(exchanges, list):
                errors.append(f"Model '{model_name}': 'exchanges' must be a list")
            elif len(exchanges) == 0:
                warnings.append(f"Model '{model_name}': 'exchanges' list is empty")
            else:
                for idx, exchange in enumerate(exchanges):
                    if not isinstance(exchange, dict):
                        errors.append(f"Model '{model_name}': exchange {idx} must be a dictionary")
                        continue
                    
                    # Check required exchange fields
                    required_fields = ["substrate", "fba_flux", "Km", "Vmax"]
                    for field in required_fields:
                        if field not in exchange:
                            errors.append(f"Model '{model_name}': exchange {idx} missing required field '{field}'")
                    
                    # Validate substrate reference
                    substrate = exchange.get("substrate")
                    if substrate:
                        all_exchange_substrates.add(substrate)
                        if substrate_names and substrate not in substrate_names:
                            errors.append(
                                f"Model '{model_name}': exchange {idx} references substrate '{substrate}' "
                                f"which is not defined in the 'substrates' section"
                            )
                    
                    # Validate numeric fields
                    for field in ["Km", "Vmax"]:
                        if field in exchange:
                            try:
                                val = float(exchange[field])
                                if val < 0:
                                    errors.append(f"Model '{model_name}': exchange {idx} '{field}' must be >= 0")
                            except (ValueError, TypeError):
                                errors.append(f"Model '{model_name}': exchange {idx} '{field}' must be a number")
        
        # Validate growth_model if present
        if "growth_model" in model_config:
            growth = model_config["growth_model"]
            if not isinstance(growth, dict):
                errors.append(f"Model '{model_name}': 'growth_model' must be a dictionary")
            else:
                numeric_fields = ["cell_density", "reference_volume", "max_growth_rate", "nuclear_volume"]
                for field in numeric_fields:
                    if field in growth:
                        try:
                            val = float(growth[field])
                            if val < 0:
                                errors.append(f"Model '{model_name}': growth_model '{field}' must be >= 0")
                        except (ValueError, TypeError):
                            errors.append(f"Model '{model_name}': growth_model '{field}' must be a number")
                
                if "objective_reaction" in growth:
                    if not isinstance(growth["objective_reaction"], str):
                        errors.append(f"Model '{model_name}': growth_model 'objective_reaction' must be a string")
        
        # Validate settings if present
        if "settings" in model_config:
            settings = model_config["settings"]
            if not isinstance(settings, dict):
                errors.append(f"Model '{model_name}': 'settings' must be a dictionary")
            else:
                if "intracellular_dt" in settings:
                    try:
                        dt = float(settings["intracellular_dt"])
                        if dt <= 0:
                            errors.append(f"Model '{model_name}': settings 'intracellular_dt' must be > 0")
                    except (ValueError, TypeError):
                        errors.append(f"Model '{model_name}': settings 'intracellular_dt' must be a number")
        
        # Validate death_model if present
        if "death_model" in model_config:
            death = model_config["death_model"]
            if not isinstance(death, dict):
                errors.append(f"Model '{model_name}': 'death_model' must be a dictionary")
            else:
                if "enabled" in death:
                    if not isinstance(death["enabled"], bool):
                        errors.append(f"Model '{model_name}': death_model 'enabled' must be a boolean")
                if "death_type" in death:
                    if not isinstance(death["death_type"], str):
                        errors.append(f"Model '{model_name}': death_model 'death_type' must be a string")
                for field in ["death_trigger_flux", "death_flux_threshold", "death_rate_increase"]:
                    if field in death:
                        try:
                            float(death[field])
                        except (ValueError, TypeError):
                            errors.append(f"Model '{model_name}': death_model '{field}' must be a number")
    
    
    # Check for substrates in exchanges that aren't defined
    if all_exchange_substrates:
        if not substrate_names:
            warnings.append(
                f"No 'substrates' section found. The following substrates are referenced in exchanges "
                f"and will use default values: {', '.join(sorted(all_exchange_substrates))}"
            )
        else:
            undefined_substrates = all_exchange_substrates - substrate_names
            if undefined_substrates:
                errors.append(
                    f"The following substrates are referenced in exchanges but not defined in 'substrates' section: "
                    f"{', '.join(sorted(undefined_substrates))}"
                )
    
    # Report warnings
    if warnings:
        print("[WARNING] Validation warnings:")
        for warning in warnings:
            print(f"  - {warning}")
    
    # Report errors and raise if any
    if errors:
        error_msg = f"YAML configuration validation failed ({len(errors)} error(s)):\n"
        for error in errors:
            error_msg += f"  - {error}\n"
        raise ValueError(error_msg)
    
    print(f"[OK] YAML configuration validation passed")
    if substrate_names:
        print(f"[INFO] Found {len(substrate_names)} substrate(s) and {len(config['models'])} model(s)")


# ----------------------
# Main orchestration
# ----------------------
def update_config_with_dfba(template_xml_path,
                           output_xml_path,
                           config_yaml_path,
                           sbml_folder,
                           template_cell_name,
                           cell_prefix="",
                           keep_existing_cells=False):
    """
    Update PhysiCell config with dFBA models from YAML configuration.
    
    Args:
        template_xml_path: Path to PhysiCell XML template
        output_xml_path: Path to write updated XML
        config_yaml_path: Path to YAML configuration file
        sbml_folder: Folder containing SBML files (for resolving relative paths)
        template_cell_name: Name of template cell_definition to clone
        cell_prefix: Optional prefix for new cell definition names
        keep_existing_cells: If False (default), remove all existing cell_definitions before adding new ones.
                            If True, keep existing cell_definitions and add new ones.
    """
    # Parse template XML
    parser = etree.XMLParser(remove_blank_text=True)
    tree = etree.parse(template_xml_path, parser)
    root = tree.getroot()

    # Load YAML configuration
    with open(config_yaml_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Validate configuration before processing
    validate_yaml_config(config, config_yaml_path, sbml_folder)
    
    models = config["models"]
    
    # Build substrate configuration dictionary
    substrate_configs = {}
    if "substrates" in config:
        for sub_config in config["substrates"]:
            name = sub_config.get("name")
            if name:
                substrate_configs[name] = {
                    "diffusion_coefficient": sub_config.get("diffusion_coefficient"),
                    "decay_rate": sub_config.get("decay_rate"),
                    "initial_condition": sub_config.get("initial_condition", 0.0)
                }
    

    # Find template cell_definition
    template_cell = find_template_cell(root, template_cell_name)

    # We'll insert new cell_definitions after the template inside <cell_definitions>
    cdefs_parent = root.find(".//cell_definitions")
    if cdefs_parent is None:
        raise ValueError("No <cell_definitions> parent found in template XML.")

    # Handle existing cell_definitions based on flag and determine starting ID
    # Note: We need to keep the template cell temporarily for cloning, so we'll remove it after cloning
    template_cell_name_to_remove = None
    next_id = 0
    
    if not keep_existing_cells:
        # Store template cell name to remove it later (after we've cloned it for all models)
        template_cell_name_to_remove = template_cell.get("name")
        # Remove all other existing cell_definitions (keep template for now so we can clone it)
        existing_cells = cdefs_parent.findall("cell_definition")
        cells_to_remove = [cell for cell in existing_cells if cell.get("name") != template_cell_name_to_remove]
        for cell in cells_to_remove:
            cdefs_parent.remove(cell)
        print(f"[INFO] Removed {len(cells_to_remove)} existing cell_definition(s) from template (keeping template for cloning)")
        # Start IDs from 0 when removing existing cells
        next_id = 0
    else:
        # Count existing cell_definitions to determine next ID
        existing_cells = cdefs_parent.findall("cell_definition")
        next_id = len(existing_cells)
        print(f"[INFO] Keeping existing cell_definitions (found {len(existing_cells)} cell_definition(s)), next ID will be {next_id}")

    created_cell_names = []
    global_substrates = set()

    # For each model: clone template and customize
    for model_name, model_config in models.items():
        # Resolve SBML path
        sbml_path_rel = model_config.get("sbml_path")
        if not sbml_path_rel:
            raise ValueError(f"Model '{model_name}' missing 'sbml_path' in configuration")
        
        # If relative path, resolve against sbml_folder
        if os.path.isabs(sbml_path_rel):
            sbml_path = sbml_path_rel
        else:
            sbml_path = os.path.join(sbml_folder, sbml_path_rel)
        
        # Check if SBML file exists
        if not os.path.exists(sbml_path):
            raise FileNotFoundError(f"SBML file not found: {sbml_path}")
        
        # Use absolute path for XML
        sbml_path_abs = os.path.abspath(sbml_path)
        
        new_name = f"{cell_prefix}{model_name}"

        template_info = f"from template '{template_cell_name}'" if template_cell_name else "from first template cell_definition"
        print(f"[INFO] Creating cell_definition '{new_name}' {template_info} with SBML: {sbml_path_abs}")

        # deep copy the template element
        new_cell = copy.deepcopy(template_cell)
        # assign new name
        new_cell.set("name", new_name)
        # assign incremental ID
        new_cell.set("ID", str(next_id))
        print(f"[INFO] Assigned ID={next_id} to cell_definition '{new_name}'")
        next_id += 1

        # Update volume block from YAML growth_model (reference_volume -> total, nuclear_volume -> nuclear)
        growth_cfg = model_config.get("growth_model", {})
        total_volume = growth_cfg.get("reference_volume")
        nuclear_volume = growth_cfg.get("nuclear_volume", 0.0)
        if total_volume is not None:
            set_cell_volume(new_cell, total_volume=total_volume, nuclear_volume=nuclear_volume)

        # Add dfba intracellular according to model configuration
        add_intracellular_dfba(new_cell, sbml_path_abs, model_config)

        # Collect all substrates from exchanges
        exchanges = model_config.get("exchanges", [])
        for exchange in exchanges:
            global_substrates.add(exchange["substrate"])

        # Ensure secretion entries for all substrates in full parameter table (global list)
        # We'll gather global list later, but for now accumulate created cell names and add later
        cdefs_parent.append(new_cell)
        created_cell_names.append(new_name)

    # If we're removing existing cells, also remove the template cell now (after cloning)
    if template_cell_name_to_remove:
        template_cell_to_remove = cdefs_parent.find(f".//cell_definition[@name='{template_cell_name_to_remove}']")
        if template_cell_to_remove is not None:
            cdefs_parent.remove(template_cell_to_remove)
            print(f"[INFO] Removed template cell_definition '{template_cell_name_to_remove}' after cloning")

    # Now we have appended new cell_definitions; compute full list of cell names (existing + created)
    all_cell_defs = [cd.get("name") for cd in cdefs_parent.findall("cell_definition")]
    all_cell_names = list(dict.fromkeys(all_cell_defs))  # preserve order, deduplicate

    # Ensure microenvironment contains all substrates
    microenv = root.find(".//microenvironment_setup")
    if microenv is None:
        raise ValueError("Template missing <microenvironment_setup> section.")
    existing_vars = {v.get("name") for v in microenv.findall("variable")}
    global_substrates_sorted = sorted(global_substrates)
    for s in global_substrates_sorted:
        if s not in existing_vars:
            # Get substrate-specific configuration or use hardcoded defaults
            sub_config = substrate_configs.get(s, {})
            diff = sub_config.get("diffusion_coefficient", 50000.0)
            decay = sub_config.get("decay_rate", 0.0)
            init = sub_config.get("initial_condition", 0.0)
            
            print(f"[INFO] Adding microenvironment variable for substrate '{s}' "
                  f"(diffusion={diff}, decay={decay}, init={init})")
            microenv.insert(0, make_microenv_variable(s, diff=diff, decay=decay, init=init))
            existing_vars.add(s)

    # For each cell_def (both original template and new ones) ensure secretion and chemotactic sensitivities and cell interactions
    for cell in cdefs_parent.findall("cell_definition"):
        name = cell.get("name")
        # ensure secretion contains all substrates
        for s in global_substrates_sorted:
            ensure_secretion_has_substrate(cell, s)
        # ensure chemotactic sensitivities contain all substrates
        for s in global_substrates_sorted:
            ensure_chemotactic_sensitivities(cell, s)
        # ensure cell_interactions entries exist for all cell types
        ensure_cell_interactions_for_all(cell, all_cell_names)

    # Finally write output
    xml_pretty_write(tree, output_xml_path)
    print(f"[OK] Wrote updated config to: {output_xml_path}")


# ----------------------
# CLI
# ----------------------
def main():
    parser = argparse.ArgumentParser(description="Inject dfBA blocks into PhysiCell config using a template cell_definition.")
    parser.add_argument("--template", "-t", required=True, help="Path to PhysiCell XML template.")
    parser.add_argument("--output", "-o", required=True, help="Path to write updated PhysiCell XML.")
    parser.add_argument("--config", "-c", required=True,
                        help="YAML configuration file with models, SBML paths, exchanges, and parameters.")
    parser.add_argument("--sbml-folder", "-s", required=True, help="Folder containing SBML files (for resolving relative paths).")
    parser.add_argument("--template-cell", default="", help="Name of the cell_definition in template to use as blueprint. If empty or not provided, uses the first cell_definition found.")
    parser.add_argument("--cell-prefix", default="", help="Optional prefix to prepend to new cell_definition names.")
    parser.add_argument("--keep-existing-cells", action="store_true", 
                        help="Keep existing cell_definitions from template. By default, all existing cell_definitions are removed before adding new ones.")
    args = parser.parse_args()

    update_config_with_dfba(
        args.template,
        args.output,
        args.config,
        args.sbml_folder,
        args.template_cell,
        cell_prefix=args.cell_prefix,
        keep_existing_cells=args.keep_existing_cells
    )


if __name__ == "__main__":
    main()

