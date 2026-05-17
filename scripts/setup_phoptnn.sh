#!/usr/bin/env bash

set -euo pipefail

INSTALL_DIR="${1:?Usage: setup_phoptnn.sh <install_dir>}"
REPO_URL="https://github.com/kuenzelab/pHoptNN.git"

mkdir -p "$(dirname "${INSTALL_DIR}")"

if [ -d "${INSTALL_DIR}/.git" ]; then
    echo "[INFO] pHoptNN repo already cloned at ${INSTALL_DIR}"
else
    echo "[INFO] Cloning pHoptNN repository ..."
    git clone "${REPO_URL}" "${INSTALL_DIR}"
fi

INTERFACE="${INSTALL_DIR}/phoptnn_interface.py"
if [ ! -f "${INTERFACE}" ]; then
    echo "[ERROR] phoptnn_interface.py not found at ${INTERFACE}" >&2
    exit 1
fi

PREDICT="${INSTALL_DIR}/EGNN/predict.py"
if [ ! -f "${PREDICT}" ]; then
    echo "[ERROR] EGNN/predict.py not found at ${PREDICT}" >&2
    exit 1
fi

W="${INSTALL_DIR}/EGNN/weight/W_6_attn.pt"
if [ -f "$W" ]; then
    echo "[INFO] Found model weights: $W"
else
    echo "[WARN] No pre-trained weights found at $W"
fi

touch "${INSTALL_DIR}/.setup_done"
echo "[INFO] pHoptNN setup complete at ${INSTALL_DIR}"
