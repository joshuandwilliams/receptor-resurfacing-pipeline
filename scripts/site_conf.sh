#!/bin/bash
#
# Load site.conf, without letting it clobber anything already set.
#
#   source scripts/site_conf.sh
#
# Precedence, highest first:
#   1. a value already in the environment, including one exported by sbatch
#      from the submitting shell
#   2. site.conf
#   3. the default each consumer applies when the value is empty
#
# That ordering is what lets a one-off run override a setting on the command
# line without editing site.conf, which matters because sbatch exports the
# submitting environment into the job.
#
# Every RECEPTOR_* name is protected, discovered with compgen rather than
# listed, so a setting added to site.conf later gets the same precedence
# without this file changing. Checked under bash 3.2, which is what macOS
# ships, and 5.x.
#
# Sourcing this when site.conf does not exist is not an error. A fresh clone
# has no site.conf, and every consumer has either a neutral default or a clear
# failure, so the scripts still run.

# Resolve the repository root from this file's location, so the launchers do
# not each need their own copy of that logic.
_SITE_CONF_SELF="${BASH_SOURCE[0]:-$0}"
_SITE_CONF_DIR="$(cd "$(dirname "${_SITE_CONF_SELF}")/.." && pwd)"
_SITE_CONF_FILE="${_SITE_CONF_DIR}/site.conf"

if [ -f "${_SITE_CONF_FILE}" ]; then
    # Snapshot every RECEPTOR_* name that is already set, with its value, then
    # restore them after sourcing. Sourcing is the only way to allow a
    # site.conf value to reference another, and it is also what would
    # otherwise let the file override the environment.
    _SITE_CONF_SAVED=""
    for _n in $(compgen -v RECEPTOR_ 2>/dev/null); do
        # Only non-empty values win. An exported-but-empty variable is treated
        # as unset, so `RECEPTOR_PARTITION_GPU= sbatch ...` does not silently
        # defeat site.conf.
        if [ -n "${!_n}" ]; then
            _SITE_CONF_SAVED="${_SITE_CONF_SAVED} ${_n}=$(printf '%q' "${!_n}")"
        fi
    done

    # shellcheck source=/dev/null
    . "${_SITE_CONF_FILE}"

    if [ -n "${_SITE_CONF_SAVED}" ]; then
        eval "export ${_SITE_CONF_SAVED}"
    fi
    unset _SITE_CONF_SAVED _n
fi

# Export everything site.conf set, so a child process and Nextflow's
# System.getenv() both see it. Without this a value assigned in site.conf is a
# shell variable only, and nextflow.config reads nothing.
for _n in $(compgen -v RECEPTOR_ 2>/dev/null); do
    export "${_n?}"
done
unset _n

# The repository root, if nothing else established it. SLURM_SUBMIT_DIR is the
# directory sbatch was called from, which is the usual case.
export RECEPTOR_PIPELINE_DIR="${RECEPTOR_PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-${_SITE_CONF_DIR}}}"

unset _SITE_CONF_SELF _SITE_CONF_DIR _SITE_CONF_FILE
