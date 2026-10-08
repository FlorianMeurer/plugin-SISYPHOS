from __future__ import annotations

import argparse
import json
import math
import re
import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


FLOAT_WITH_ESD_RE = re.compile(
    r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)(?:\((\d+)\))?$"
)


@dataclass
class ParsedNumber:
    value: float
    esd: float
    raw: str


def _split_mantissa_and_exp(number_str: str) -> Tuple[str, int]:
    """Split a numeric string into mantissa and base-10 exponent."""
    lowered = number_str.lower()
    if "e" in lowered:
        mantissa, exp = lowered.split("e", 1)
        return mantissa, int(exp)
    return number_str, 0


def parse_cif_number_with_esd(token: str) -> ParsedNumber:
    """Parse CIF values like 0.023(4), 1.2e-2(3), ., ? into value + esd."""
    cleaned = token.strip()
    if cleaned in {".", "?", ""}:
        return ParsedNumber(value=math.nan, esd=math.nan, raw=token)

    match = FLOAT_WITH_ESD_RE.match(cleaned)
    if not match:
        raise ValueError(f"Unsupported numeric token format: '{token}'")

    value_str = match.group(1)
    esd_digits_str = match.group(2)

    value = float(value_str)

    if esd_digits_str is None:
        return ParsedNumber(value=value, esd=math.nan, raw=token)

    mantissa, exponent = _split_mantissa_and_exp(value_str)
    decimals = len(mantissa.split(".", 1)[1]) if "." in mantissa else 0
    esd_digits = int(esd_digits_str)
    esd = esd_digits * (10 ** (exponent - decimals))

    return ParsedNumber(value=value, esd=esd, raw=token)


def _tokenize_cif_row(line: str) -> List[str]:
    """Tokenize one CIF row with quoted-value support."""
    stripped = line.strip()
    if not stripped:
        return []
    return shlex.split(stripped, comments=False, posix=True)


def _collect_loop(lines: List[str], start_idx: int) -> Tuple[List[str], List[List[str]], int]:
    """Collect a CIF loop_ block and return headers, data rows, and next index."""
    idx = start_idx + 1
    headers: List[str] = []

    while idx < len(lines):
        current = lines[idx].strip()
        if current.startswith("_"):
            headers.append(current)
            idx += 1
            continue
        break

    rows: List[List[str]] = []
    while idx < len(lines):
        raw_line = lines[idx]
        stripped = raw_line.strip()

        if not stripped or stripped.startswith("#"):
            idx += 1
            continue

        # End of current loop data.
        if stripped.startswith(("_", "loop_", "data_", "save_")):
            break

        # Multi-line text values starting with ';' are not expected for atom rows.
        if stripped.startswith(";"):
            idx += 1
            while idx < len(lines) and not lines[idx].startswith(";"):
                idx += 1
            if idx < len(lines):
                idx += 1
            continue

        tokens = _tokenize_cif_row(raw_line)
        if tokens:
            rows.append(tokens)
        idx += 1

    return headers, rows, idx


def _header_index_map(headers: List[str]) -> Dict[str, int]:
    return {header.lower(): i for i, header in enumerate(headers)}


def extract_u_values_from_cif(cif_path: str) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Extract Uiso and anisotropic U components from CIF loops.

    Returns:
        {
          "u_iso": {
            "C1": {"value": 0.0123, "esd": 0.0004, "raw": "0.0123(4)"},
            ...
          },
          "u_aniso": {
            "C1": {
              "u11": {"value": ..., "esd": ..., "raw": ...},
              "u22": {...},
              "u33": {...},
              "u12": {...},
              "u13": {...},
              "u23": {...}
            },
            ...
          }
        }
    """
    path = Path(cif_path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    result: Dict[str, Dict[str, Dict[str, float]]] = {
        "u_iso": {},
        "u_aniso": {},
    }

    idx = 0
    while idx < len(lines):
        line = lines[idx].strip()

        if line != "loop_":
            idx += 1
            continue

        headers, rows, next_idx = _collect_loop(lines, idx)
        header_map = _header_index_map(headers)

        # Atom-site loop with Uiso/Ueq.
        atom_label_key = "_atom_site_label"
        atom_u_iso_key = "_atom_site_u_iso_or_equiv"

        if atom_label_key in header_map and atom_u_iso_key in header_map:
            i_label = header_map[atom_label_key]
            i_u_iso = header_map[atom_u_iso_key]

            for row in rows:
                if len(row) <= max(i_label, i_u_iso):
                    continue

                label = row[i_label]
                parsed = parse_cif_number_with_esd(row[i_u_iso])
                result["u_iso"][label] = {
                    "value": parsed.value,
                    "esd": parsed.esd,
                    "raw": parsed.raw,
                }

        # Anisotropic displacement loop.
        aniso_label_key = "_atom_site_aniso_label"
        aniso_keys = {
            "u11": "_atom_site_aniso_u_11",
            "u22": "_atom_site_aniso_u_22",
            "u33": "_atom_site_aniso_u_33",
            "u12": "_atom_site_aniso_u_12",
            "u13": "_atom_site_aniso_u_13",
            "u23": "_atom_site_aniso_u_23",
        }

        if aniso_label_key in header_map and all(k in header_map for k in aniso_keys.values()):
            i_label = header_map[aniso_label_key]
            i_u = {name: header_map[cif_key] for name, cif_key in aniso_keys.items()}

            for row in rows:
                required_max = max([i_label] + list(i_u.values()))
                if len(row) <= required_max:
                    continue

                label = row[i_label]
                components: Dict[str, Dict[str, float]] = {}

                for comp_name, comp_idx in i_u.items():
                    parsed = parse_cif_number_with_esd(row[comp_idx])
                    components[comp_name] = {
                        "value": parsed.value,
                        "esd": parsed.esd,
                        "raw": parsed.raw,
                    }

                result["u_aniso"][label] = components

        idx = next_idx

    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract Uiso and anisotropic U values (with ESDs) from a CIF file."
    )
    parser.add_argument("cif", help="Path to the CIF file")
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="JSON indentation (default: 2)",
    )
    args = parser.parse_args()

    extracted = extract_u_values_from_cif(args.cif)
    print(json.dumps(extracted, indent=args.indent, allow_nan=True))


if __name__ == "__main__":
    main()
