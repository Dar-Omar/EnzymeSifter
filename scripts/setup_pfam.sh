#!/usr/bin/env bash

set -euo pipefail

INSTALL_DIR="${1:?Usage: setup_pfam.sh <install_dir>}"
URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
GZ_FILE="${INSTALL_DIR}/Pfam-A.hmm.gz"
HMM_FILE="${INSTALL_DIR}/Pfam-A.hmm"
MAX_RETRIES=2

mkdir -p "${INSTALL_DIR}"

download_ok=false
for attempt in $(seq 0 "${MAX_RETRIES}"); do
    if [ -f "${GZ_FILE}" ]; then
        echo "[INFO] Verifying existing Pfam-A.hmm.gz (attempt ${attempt})..."
        if gzip -t "${GZ_FILE}" 2>/dev/null; then
            echo "[INFO] Archive integrity OK."
            download_ok=true
            break
        else
            echo "[WARN] Archive is corrupt. Deleting and re-downloading..."
            rm -f "${GZ_FILE}"
        fi
    fi

    if [ "${attempt}" -gt 0 ]; then
        echo "[INFO] Retry ${attempt}/${MAX_RETRIES}..."
    fi
    echo "[INFO] Downloading Pfam-A.hmm.gz (~1.5 GB) from EBI..."

    if command -v wget &>/dev/null; then
        wget --no-verbose --tries=3 --timeout=120 \
             -O "${GZ_FILE}" "${URL}"
    elif command -v curl &>/dev/null; then
        curl -L --retry 3 --connect-timeout 120 \
             -o "${GZ_FILE}" "${URL}"
    else
        echo "[ERROR] Neither wget nor curl found." >&2
        exit 1
    fi
done

if [ "${download_ok}" = false ]; then
    if gzip -t "${GZ_FILE}" 2>/dev/null; then
        echo "[INFO] Archive integrity OK after download."
    else
        echo "[ERROR] Archive still corrupt after ${MAX_RETRIES} retries." >&2
        exit 1
    fi
fi

if [ ! -f "${HMM_FILE}" ]; then
    echo "[INFO] Decompressing Pfam-A.hmm..."
    gunzip -k "${GZ_FILE}"
else
    echo "[INFO] Pfam-A.hmm already decompressed."
fi

SSI_FILE="${HMM_FILE}.ssi"
STRIPPED_SENTINEL="${INSTALL_DIR}/.acc_stripped"

if [ ! -f "${SSI_FILE}" ]; then
    if ! command -v hmmfetch &>/dev/null; then
        echo "[ERROR] hmmfetch not found on PATH. Is the HMMER conda env active?" >&2
        exit 1
    fi

    if [ ! -f "${STRIPPED_SENTINEL}" ]; then
        echo "[INFO] Stripping ACC version suffixes (PFXXXXX.N -> PFXXXXX) "
        echo "       so hmmfetch can resolve bare accessions..."
        sed -i 's/^\(ACC[[:space:]]\+PF[0-9]\{5\}\)\.[0-9]\+/\1/' "${HMM_FILE}"
        touch "${STRIPPED_SENTINEL}"
    fi

    echo "[INFO] Building hmmfetch index (Pfam-A.hmm.ssi)..."
    hmmfetch --index "${HMM_FILE}"
else
    echo "[INFO] hmmfetch index already present."
fi

touch "${INSTALL_DIR}/.setup_done"
echo "[INFO] Pfam setup complete at ${INSTALL_DIR}"
