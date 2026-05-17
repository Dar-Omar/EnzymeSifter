#!/usr/bin/env bash

set -euo pipefail

INSTALL_DIR="${1:?Usage: setup_netsolp.sh <install_dir>}"
URL="https://services.healthtech.dtu.dk/services/NetSolP-1.0/netsolp-1.0.ALL.tar.gz"
TARBALL="${INSTALL_DIR}/netsolp-1.0.ALL.tar.gz"
MAX_RETRIES=2

mkdir -p "${INSTALL_DIR}"

download_ok=false
for attempt in $(seq 0 "${MAX_RETRIES}"); do
    if [ -f "${TARBALL}" ]; then
        echo "[INFO] Verifying existing tarball (attempt ${attempt})..."
        if gzip -t "${TARBALL}" 2>/dev/null; then
            echo "[INFO] Tarball integrity OK."
            download_ok=true
            break
        else
            echo "[WARN] Tarball is corrupt. Deleting and re-downloading..."
            rm -f "${TARBALL}"
        fi
    fi

    if [ "${attempt}" -gt 0 ]; then
        echo "[INFO] Retry ${attempt}/${MAX_RETRIES}..."
    fi
    echo "[INFO] Downloading NetSolP-1.0 (~5.6 GB) from DTU ..."

    if command -v wget &>/dev/null; then
        wget --no-verbose --tries=3 --timeout=120 \
             -O "${TARBALL}" "${URL}"
    elif command -v curl &>/dev/null; then
        curl -L --retry 3 --connect-timeout 120 \
             -o "${TARBALL}" "${URL}"
    else
        echo "[ERROR] Neither wget nor curl found." >&2
        exit 1
    fi
done

if [ "${download_ok}" = false ]; then
    if gzip -t "${TARBALL}" 2>/dev/null; then
        echo "[INFO] Tarball integrity OK after download."
    else
        echo "[ERROR] Tarball still corrupt after ${MAX_RETRIES} retries." >&2
        exit 1
    fi
fi

echo "[INFO] Extracting PredictionServer code ..."
tar -xzf "${TARBALL}" \
    -C "${INSTALL_DIR}" \
    --strip-components=1 \
    "Data_and_Code/PredictionServer"

echo "[INFO] Extracting ONNX models into PredictionServer/models/ ..."
tar -xzf "${TARBALL}" \
    -C "${INSTALL_DIR}/PredictionServer" \
    "models"

if tar -tzf "${TARBALL}" "data.py" &>/dev/null; then
    echo "[INFO] Extracting top-level data.py ..."
    tar -xzf "${TARBALL}" -C "${INSTALL_DIR}/PredictionServer" "data.py"
fi

PREDICT_PY="${INSTALL_DIR}/PredictionServer/predict.py"
if [ ! -f "${PREDICT_PY}" ]; then
    echo "[ERROR] predict.py not found at ${PREDICT_PY}" >&2
    exit 1
fi

ONNX_COUNT=$(find "${INSTALL_DIR}/PredictionServer/models" -name "*.onnx" 2>/dev/null | wc -l)
echo "[INFO] NetSolP setup complete.  predict.py → ${PREDICT_PY}"
echo "[INFO] Found ${ONNX_COUNT} ONNX model file(s)."

if [ "${ONNX_COUNT}" -eq 0 ]; then
    echo "[ERROR] No ONNX models found — prediction will fail." >&2
    exit 1
fi

touch "${INSTALL_DIR}/PredictionServer/.setup_done"


# Clean up duplicated files
rm -f "${INSTALL_DIR}/"*_quantized.onnx
rm -f "${INSTALL_DIR}/ESM"*"_alphabet.pkl"
rm -f "${INSTALL_DIR}/PredictionServer/"*_quantized.onnx
rm -f "${INSTALL_DIR}/PredictionServer/ESM"*"_alphabet.pkl"
