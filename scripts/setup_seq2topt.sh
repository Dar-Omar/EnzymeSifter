#!/usr/bin/env bash

set -euo pipefail

INSTALL_DIR="${1:?Usage: setup_seq2topt.sh <install_dir>}"
REPO_URL="https://github.com/SizheQiu/Seq2Topt.git"
RELEASE_URL="https://github.com/SizheQiu/Seq2Topt/releases/download/v1.0.0"

mkdir -p "$(dirname "${INSTALL_DIR}")"

if [ -d "${INSTALL_DIR}/.git" ]; then
    echo "[INFO] Seq2Topt repo already cloned at ${INSTALL_DIR}"
else
    echo "[INFO] Cloning Seq2Topt repository ..."
    git clone "${REPO_URL}" "${INSTALL_DIR}"
fi

for F in "code/seq2topt.py" "code/model.py"; do
    if [ ! -f "${INSTALL_DIR}/${F}" ]; then
        echo "[ERROR] ${F} not found at ${INSTALL_DIR}/${F}" >&2
        exit 1
    fi
done

WEIGHTS_DIR="${INSTALL_DIR}/weights"
mkdir -p "${WEIGHTS_DIR}"

for W in model_topt_window.3_r2.0.57.pth model_tm_window.3_r2.0.76.pth; do
    DEST="${WEIGHTS_DIR}/${W}"
    if [ -f "${DEST}" ]; then
        echo "[INFO] Weight file already exists: ${DEST}"
    else
        echo "[INFO] Downloading ${W} ..."
        if command -v wget &>/dev/null; then
            wget --no-verbose -O "${DEST}" "${RELEASE_URL}/${W}"
        elif command -v curl &>/dev/null; then
            curl -L -o "${DEST}" "${RELEASE_URL}/${W}"
        else
            echo "[ERROR] Neither wget nor curl found." >&2
            exit 1
        fi

        if [ ! -s "${DEST}" ]; then
            echo "[ERROR] Downloaded weight file is empty: ${DEST}" >&2
            rm -f "${DEST}"
            exit 1
        fi
        echo "[INFO] Downloaded: ${DEST}"
    fi
done

touch "${INSTALL_DIR}/.setup_done"
echo "[INFO] Seq2Topt setup complete at ${INSTALL_DIR}"
