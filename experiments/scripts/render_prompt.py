#!/usr/bin/env python3
"""
render_prompt.py — Render the literature-sweep prompt template for a campaign.

Reads a campaign config (YAML) and a prompt template (Markdown with placeholder
syntax), substitutes the config values, and writes the rendered prompt out.

Template syntax:
    {{FIELD}}                       Required substitution. Errors if FIELD missing.
    {{?BLOCK_NAME}}...{{/BLOCK_NAME}}   Optional block. Included only if BLOCK_NAME
                                    is present and non-empty in the config. The
                                    block body may itself contain {{FIELD}}
                                    substitutions which are resolved when included.

Config schema (YAML):
    receptor:
        name: <str>                 REQUIRED
        uniprot_id: <str>           optional
        domain_boundaries: <str>    optional (e.g. "186-263")
    effector:
        name: <str>                 REQUIRED
        aliases: [<str>, ...]       optional
    campaign_slug: <str>            optional (used in output filenames)
    prior_knowledge: [<str>, ...]   optional (free-text bullets)

Usage:
    python render_prompt.py CONFIG.yml TEMPLATE.md.template [-o OUTPUT.md]

Exits non-zero on missing required fields or template/config parse errors.
"""

import argparse
import re
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    print(
        "ERROR: PyYAML is required. Install with: pip install pyyaml",
        file=sys.stderr,
    )
    sys.exit(2)


# ── Required-field map ────────────────────────────────────────────────────────
# Maps the placeholder name (as it appears in {{FIELD}}) to the dotted path
# in the YAML config. Anything in this map MUST be present in the config.

REQUIRED_FIELDS = {
    "RECEPTOR_NAME": "receptor.name",
    "EFFECTOR_NAME": "effector.name",
}

# Optional-field map: placeholder name → dotted path. Substituted if present.
OPTIONAL_FIELDS = {
    "RECEPTOR_UNIPROT": "receptor.uniprot_id",
    "RECEPTOR_DOMAIN_BOUNDARIES": "receptor.domain_boundaries",
    "EFFECTOR_ALIASES": "effector.aliases",
    "CAMPAIGN_SLUG": "campaign_slug",
    "PRIOR_KNOWLEDGE": "prior_knowledge",
}

# Optional-block map: block name → dotted path. The block is included only if
# the path resolves to a non-empty value in the config.
OPTIONAL_BLOCKS = {
    "HAS_RECEPTOR_UNIPROT": "receptor.uniprot_id",
    "HAS_RECEPTOR_DOMAIN_BOUNDARIES": "receptor.domain_boundaries",
    "HAS_EFFECTOR_ALIASES": "effector.aliases",
    "HAS_PRIOR_KNOWLEDGE": "prior_knowledge",
}


# ── Helpers ───────────────────────────────────────────────────────────────────


def resolve_path(config: dict, dotted_path: str):
    """Resolve a dotted path against a nested config dict. Returns None if any
    intermediate key is missing."""
    node = config
    for key in dotted_path.split("."):
        if not isinstance(node, dict) or key not in node:
            return None
        node = node[key]
    return node


def is_nonempty(value) -> bool:
    """A field is 'present' if it's not None, not empty string, not empty list."""
    if value is None:
        return False
    if isinstance(value, (str, list, dict)) and len(value) == 0:
        return False
    return True


def format_field_value(value) -> str:
    """Format a field value for substitution into the prompt body.

    Strings: as-is.
    Lists: rendered as a Markdown bullet list, one item per line, with leading
           newline so the bullets sit on their own lines under whatever heading
           introduced them in the template.
    Other: str().
    """
    if isinstance(value, str):
        return value
    if isinstance(value, list):
        if not value:
            return ""
        return "\n" + "\n".join(f"- {item}" for item in value)
    return str(value)


# ── Block expansion ───────────────────────────────────────────────────────────

# Match {{?BLOCK_NAME}}...{{/BLOCK_NAME}} (non-greedy, multiline). The pattern
# also absorbs the newline immediately after the opening tag and the newline
# immediately before the closing tag if they exist, so block markers can sit
# on their own lines in the template without producing extra blank lines in
# the rendered output.
BLOCK_PATTERN = re.compile(
    r"\{\{\?(?P<name>[A-Z_]+)\}\}\n?(?P<body>.*?)\n?\{\{/(?P=name)\}\}",
    re.DOTALL,
)

# Match {{FIELD}}.
FIELD_PATTERN = re.compile(r"\{\{(?P<name>[A-Z_]+)\}\}")


def expand_blocks(template: str, config: dict) -> str:
    """Process {{?BLOCK_NAME}}...{{/BLOCK_NAME}} blocks. Includes the body if the
    associated config path resolves to a non-empty value, otherwise removes it."""

    def replace_block(match: re.Match) -> str:
        name = match.group("name")
        body = match.group("body")
        if name not in OPTIONAL_BLOCKS:
            print(
                f"ERROR: Unknown optional block: {{{{?{name}}}}}",
                file=sys.stderr,
            )
            sys.exit(1)
        path = OPTIONAL_BLOCKS[name]
        value = resolve_path(config, path)
        if is_nonempty(value):
            return body
        # Strip the block AND a trailing newline if present, so removed blocks
        # don't leave blank lines littering the output.
        return ""

    # Loop until no more blocks (handles potential nesting, though we don't
    # currently use it).
    while True:
        new_template, n = BLOCK_PATTERN.subn(replace_block, template)
        if n == 0:
            return new_template
        template = new_template


def expand_fields(template: str, config: dict) -> str:
    """Substitute {{FIELD}} placeholders. Errors on missing required fields;
    silently empties missing optional fields (which should already have been
    handled by block-expansion, but we're defensive)."""

    def replace_field(match: re.Match) -> str:
        name = match.group("name")
        if name in REQUIRED_FIELDS:
            value = resolve_path(config, REQUIRED_FIELDS[name])
            if not is_nonempty(value):
                print(
                    f"ERROR: Required field missing or empty: {REQUIRED_FIELDS[name]} "
                    f"(placeholder {{{{{name}}}}})",
                    file=sys.stderr,
                )
                sys.exit(1)
            return format_field_value(value)
        if name in OPTIONAL_FIELDS:
            value = resolve_path(config, OPTIONAL_FIELDS[name])
            return format_field_value(value) if is_nonempty(value) else ""
        print(
            f"ERROR: Unknown placeholder: {{{{{name}}}}}",
            file=sys.stderr,
        )
        sys.exit(1)

    return FIELD_PATTERN.sub(replace_field, template)


def tidy_blank_lines(text: str) -> str:
    """Collapse runs of 3+ blank lines (which can result from removed optional
    blocks) into 2. Cosmetic only."""
    return re.sub(r"\n{3,}", "\n\n", text)


# ── Main ──────────────────────────────────────────────────────────────────────


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Render the literature-sweep prompt template for a campaign, "
            "substituting values from a YAML config."
        ),
    )
    parser.add_argument("config", type=Path, help="Campaign config YAML file")
    parser.add_argument("template", type=Path, help="Prompt template (Markdown)")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output file (default: stdout)",
    )
    args = parser.parse_args()

    if not args.config.exists():
        print(f"ERROR: Config file not found: {args.config}", file=sys.stderr)
        return 1
    if not args.template.exists():
        print(f"ERROR: Template file not found: {args.template}", file=sys.stderr)
        return 1

    try:
        config = yaml.safe_load(args.config.read_text())
    except yaml.YAMLError as e:
        print(f"ERROR: Failed to parse YAML config: {e}", file=sys.stderr)
        return 1
    if not isinstance(config, dict):
        print(
            f"ERROR: Config must be a YAML mapping at the top level, got "
            f"{type(config).__name__}",
            file=sys.stderr,
        )
        return 1

    template = args.template.read_text()

    # Order matters: blocks first, fields second. Block bodies may contain
    # field placeholders that resolve correctly only after blocks are expanded.
    rendered = expand_blocks(template, config)
    rendered = expand_fields(rendered, config)
    rendered = tidy_blank_lines(rendered)

    if args.output:
        args.output.write_text(rendered)
        print(f"Rendered prompt written to {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(rendered)

    return 0


if __name__ == "__main__":
    sys.exit(main())