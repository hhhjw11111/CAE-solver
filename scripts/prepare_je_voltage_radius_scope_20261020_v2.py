"""Prepare the expanded 0--8 V LQ/NT and radius matrix; never launch a solver.

The previous nine-case protocol and original candidate inputs are retained.
Cold starts explicitly have zero restart load; final target heat remains one. New inputs are explicitly
candidate-grid, unsolved preparations, not mesh or production certificates.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import prepare_je_paper_input_20261018 as old
from audit_je_paper_reuse_20261018 import validate_acceptance

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT/"runs/JE_RECONSTRUCTION_20261018/paper_scope"
CONDITIONS = {(3, model, i/2) for model in ("LQ","NT") for i in range(17)}
CONDITIONS |= {(radius,"NT",voltage) for radius in (2,4) for voltage in (6.,8.)}


def make_input(template, *, radius_um, model, voltage, mesh_path, seed_path=None, seed_voltage=None):
    text = template.replace("\r\n","\n")
    if hashlib.sha256(text.encode()).hexdigest() != old.TEMPLATE_TEXT_SHA256:
        raise ValueError("Frozen JE template changed")
    if (radius_um,model,voltage) not in CONDITIONS:
        raise ValueError("Outside the expanded 38-case scope")
    if seed_path is None and seed_voltage is not None:
        raise ValueError("Seed voltage without a seed")
    if seed_path is not None:
        if model != "NT" or voltage == 0 or seed_voltage is None or not math.isfinite(seed_voltage) or not 0 <= seed_voltage <= voltage:
            raise ValueError("Only a verified same-model NT seed at 0<=Vs<=V is allowed")
    changes = dict(
        case_id=f"'JE_PAPER_R{radius_um:02d}_{model}_V{round(voltage*10):03d}'",
        mesh_file=old._quoted_path(mesh_path), material_profile=f"'ZNO_F1_{model}'",
        hole_radius=f"{radius_um*1e-6:.17g}", target_voltage=f"{-voltage:.17g}",
        target_voltage_factor="1.0", target_joule_factor="1.0",
        heat_source_current_work=".true.", transport_supg=".true." if model=="NT" else ".false.",
        restart_nodes_file="''" if seed_path is None else old._quoted_path(seed_path),
        restart_voltage_factor="0.0" if seed_path is None else f"{seed_voltage/voltage:.17g}",
        restart_joule_factor="0.0" if seed_path is None else "1.0")
    before = text
    for key,value in changes.items():
        text=old._replace(text,key,value)
    masked_before, masked_after = before, text
    for key in changes:
        masked_before=old._replace(masked_before,key,"UNCHANGED_COMPARISON")
        masked_after=old._replace(masked_after,key,"UNCHANGED_COMPARISON")
    if masked_before != masked_after:
        raise ValueError("Unexpected parameter change")
    for literal in ("residual_tolerance=1e-9,update_tolerance=1e-9",
                    "field_scales=1e-10,1e-10,1e-10,8.0,0.1,0.1,35.0",
                    "arclength_enabled=.false.", "heat_source_current_work=.true."):
        if literal not in text:
            raise ValueError("Original numerical/model setting changed")
    text.encode("ascii")
    return text


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    validate_acceptance(PAPER/"verification_reuse/acceptance.json")
    if old.sha(old.TEMPLATE) != old.TEMPLATE_SHA256:
        raise ValueError("Original template identity changed")
    template=old.TEMPLATE.read_text(encoding="utf-8-sig")
    geometry=PAPER/"production_preparation/geometry_validation_01/r3_A48_R22_Z18"
    mesh=geometry/"H1.mesh"
    old.validate_mesh(mesh,geometry/"manifest.json",3)
    native_mesh=old.native_path(mesh)
    # All combinations and negative controls are input tests, not native solves.
    checks=[]
    for radius,model,voltage in sorted(CONDITIONS):
        text=make_input(template,radius_um=radius,model=model,voltage=voltage,mesh_path="C:/input_test/H1.mesh")
        assert f"material_profile='ZNO_F1_{model}'" in text
        assert ("transport_supg=.true." if model=="NT" else "transport_supg=.false.") in text
        checks.append(dict(radius_um=radius,model=model,voltage_V=voltage,input_logic_passed=True))
    for invalid in [(3,"N0",8),(3,"NT",8.5),(3,"LQ",.25),(2,"LQ",6),(4,"NT",7)]:
        try:
            make_input(template,radius_um=invalid[0],model=invalid[1],voltage=invalid[2],mesh_path="C:/input_test/H1.mesh")
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid scope accepted")
    try:
        make_input(template,radius_um=3,model="NT",voltage=0,mesh_path="C:/input_test/H1.mesh",
                   seed_path="C:/input_test/nodes.csv",seed_voltage=0)
    except ValueError:
        pass
    else:
        raise AssertionError("Zero-voltage seed division was not rejected")
    args.output.mkdir(parents=True,exist_ok=False)
    records=[]
    for radius,model,voltage in sorted(CONDITIONS):
        case_id=f"R{radius:02d}_{model}_V{round(voltage*10):03d}"
        row=dict(case_id=case_id,radius_um=radius,model=model,voltage_V=voltage,
            heat_source="SIGNED_E_DOT_MODEL_CURRENT_R0",status="AWAITING_MODEL_SPECIFIC_GRID_AND_SEED",
            prepared_input="",grid_selected=False,production_acceptance=False)
        if radius==3 and (model=="LQ" or voltage==0):
            folder=args.output/"candidate_inputs"/case_id
            folder.mkdir(parents=True,exist_ok=False)
            inp=folder/"input.nml"
            inp.write_text(make_input(template,radius_um=radius,model=model,voltage=voltage,
                                     mesh_path=native_mesh),encoding="ascii")
            row.update(status="INPUT_ONLY_CANDIDATE_GRID_NOT_SOLVED",prepared_input=str(inp.resolve()))
            old.save(folder/"preparation_manifest.json",dict(**row,input_sha256=old.sha(inp),
                mesh=str(mesh),mesh_sha256=old.sha(mesh),candidate_counts=[48,22,18],
                material_profile=f"ZNO_F1_{model}",transport_supg=model=="NT",
                native_runs=0,seed="COLD_STATE_NO_CROSS_MODEL_TRANSFER",
                validation_pending=["native inspect","native solve","independent flux","model-specific spatial checks"]))
        records.append(row)
    with (args.output/"case_matrix.csv").open("x",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    old.save(args.output/"scope_contract.json",dict(
        scope="JE_PAPER_VOLTAGE_RADIUS_20261020",case_count=len(records),
        voltage_values_V=[i/2 for i in range(17)],voltage_step_assumption_V=.5,
        main_radius_um=3,voltage_models=["LQ","NT"],
        radius_comparison=dict(radii_um=[2,3,4],voltages_V=[6,8],model="NT",r3_endpoints_reused=True),
        previous_nine_case_scope_preserved=True,thresholds_changed=False,physics_changed=False,
        LQ_definition="Linearized carrier transport with frozen mobility/diffusion and three-stage one-way coupling; retains quadratic JE",
        NT_definition="Exponential carriers, temperature-dependent mobility/diffusion and fully coupled feedback",
        comparison_interpretation="Combined model-approximation difference, not isolated temperature-dependence attribution",
        NT_acceptance="Original joint mesh, SUPG, local fields and independent flux gates; all FAIL retained",
        LQ_acceptance="Own spatial and independent-flux validation; nonlinear SUPG gate is not applicable",
        radius_acceptance="Radius-specific mesh/local-field/SUPG evidence required; no inherited r3 certificate",
        zero_voltage="Native equilibrium verification; absolute/reference-scaled errors instead of division by zero power",
        main_figure_metrics=["theta_max_K","joule_W","ring_E_mean_V_m","ring_Q_mean_W_m3"],
        cloud_metrics=["signed_JE","phi","n","p","coupled_total_von_Mises"],
        geometry_comparison="Keep outer dimensions, thickness, materials, boundary loads and 0.5 um annulus width fixed",
        plotting_rules=["Only actual states; no old drift/JE splicing","No invented unsolved points","No prescribed direction of trends",
                        "Common data-driven color scales","Original field-map style retained"],
        all_results_available=False,production_acceptance=False))
    old.save(args.output/"input_logic_checks.json",dict(valid_cases=checks,valid_case_count=38,
        rejected_invalid_scope_cases=5,zero_voltage_seed_rejected=True,
        candidate_input_count=sum(bool(x["prepared_input"]) for x in records),native_runs=0))
    evidence=[Path(__file__).resolve(),Path(old.__file__).resolve(),old.TEMPLATE,
              mesh,geometry/"manifest.json",PAPER/"verification_reuse/acceptance.json"]
    old.save(args.output/"provenance.json",dict(sources={str(p):old.sha(p) for p in evidence},
        files={str(p.relative_to(args.output)):old.sha(p) for p in args.output.rglob("*") if p.is_file()},
        native_runs=0,production_acceptance=False))
    print(json.dumps(dict(cases=len(records),prepared_inputs=sum(bool(x["prepared_input"]) for x in records),
                          native_runs=0,production_acceptance=False)))


if __name__=="__main__":
    main()
