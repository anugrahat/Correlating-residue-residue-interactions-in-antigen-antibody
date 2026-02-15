"""
Residue numbering mapping for C1-87g7 WT system.

GRO sequential (from pdb2gmx):
  Chain A (antigen): GRO resid 1-194
  Chain H (heavy):   GRO resid 195-314
  Chain L (light):   GRO resid 315-421

Global numbering (from elastic net regression / setup_smd_campaign.py):
  Chain A (antigen): global 1-195    -> PDB chain A resid 1-195
  Chain L (light):   global 196-302  -> PDB chain L resid 1-107
  Chain H (heavy):   global 303-422  -> PDB chain H resid 1-120

NOTE: Global puts L before H, but GRO puts H before L!

Key mutations (global -> GRO):
  Y404 -> GRO 296, Y405 -> GRO 297, Y406 -> GRO 298  (hot spots, Chain H)
  F287 -> GRO 406  (warm spot, Chain L)
  D265 -> GRO 384  (warm spot, Chain L)
  S262 -> GRO 381  (warm spot, Chain L)
  A246 -> GRO 365  (warm spot, Chain L)
  N248 -> GRO 367  (warm spot, Chain L)
"""


def gro_to_global(gro_resid):
    """Convert GRO sequential resid to global numbering."""
    if 1 <= gro_resid <= 194:
        # Chain A (antigen): GRO 1-194 -> global 1-194
        # Note: PDB chain A starts at resid 2, but global starts at 1
        # Actually global 1-195 for antigen, but GRO only has 194 residues
        # Global for antigen: gro_resid maps to gro_resid (same)
        return gro_resid
    elif 195 <= gro_resid <= 314:
        # Chain H (heavy): GRO 195-314 -> PDB H resid 1-120 -> global 303-422
        pdb_resid = gro_resid - 194  # 1-120
        return pdb_resid + 302        # 303-422
    elif 315 <= gro_resid <= 421:
        # Chain L (light): GRO 315-421 -> PDB L resid 1-107 -> global 196-302
        pdb_resid = gro_resid - 314  # 1-107
        return pdb_resid + 195        # 196-302
    else:
        raise ValueError(f"GRO resid {gro_resid} out of range 1-421")


def global_to_gro(global_resid):
    """Convert global numbering to GRO sequential resid."""
    if 1 <= global_resid <= 195:
        # Chain A (antigen)
        return global_resid
    elif 196 <= global_resid <= 302:
        # Chain L (light): global 196-302 -> PDB L resid 1-107 -> GRO 315-421
        pdb_resid = global_resid - 195  # 1-107
        return pdb_resid + 314           # 315-421
    elif 303 <= global_resid <= 422:
        # Chain H (heavy): global 303-422 -> PDB H resid 1-120 -> GRO 195-314
        pdb_resid = global_resid - 302  # 1-120
        return pdb_resid + 194           # 195-314
    else:
        raise ValueError(f"Global resid {global_resid} out of range 1-422")


def gro_to_chain(gro_resid):
    """Convert GRO resid to (chain_letter, pdb_resid)."""
    if 1 <= gro_resid <= 194:
        return "A", gro_resid + 1  # PDB chain A starts at resid 2
    elif 195 <= gro_resid <= 314:
        return "H", gro_resid - 194
    elif 315 <= gro_resid <= 421:
        return "L", gro_resid - 314
    else:
        raise ValueError(f"GRO resid {gro_resid} out of range 1-421")


if __name__ == "__main__":
    # Verify known mutations
    test_cases = [
        (405, "Y405 (Chain H hot spot)"),
        (404, "Y404 (Chain H hot spot)"),
        (406, "Y406 (Chain H hot spot)"),
        (287, "F287 (Chain L warm spot)"),
        (265, "D265 (Chain L warm spot)"),
        (262, "S262 (Chain L warm spot)"),
        (246, "A246 (Chain L warm spot)"),
        (248, "N248 (Chain L warm spot)"),
    ]
    for global_id, name in test_cases:
        gro_id = global_to_gro(global_id)
        chain, pdb_id = gro_to_chain(gro_id)
        back = gro_to_global(gro_id)
        print(f"Global {global_id} ({name}) -> GRO {gro_id} -> Chain {chain} resid {pdb_id} -> Global {back}")
