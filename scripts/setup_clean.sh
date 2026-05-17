#!/usr/bin/env bash

set -euo pipefail

INSTALL_DIR="${1:?Usage: setup_clean.sh <install_dir>}"
REPO_URL="https://github.com/tttianhao/CLEAN.git"
ESM_REPO_URL="https://github.com/facebookresearch/esm.git"
WEIGHTS_GDRIVE_ID="1kwYd4VtzYuMvJMWXy6Vks91DSUAOcKpZ"

mkdir -p "$(dirname "${INSTALL_DIR}")"
INSTALL_DIR="$(cd "$(dirname "${INSTALL_DIR}")" && pwd)/$(basename "${INSTALL_DIR}")"

if [ -d "${INSTALL_DIR}/.git" ]; then
    echo "[INFO] CLEAN repo already cloned at ${INSTALL_DIR}"
else
    echo "[INFO] Cloning CLEAN repository ..."
    git clone "${REPO_URL}" "${INSTALL_DIR}"
fi

APP_DIR="${INSTALL_DIR}/app"
INFER_PY="${APP_DIR}/CLEAN_infer_fasta.py"
if [ ! -f "${INFER_PY}" ]; then
    echo "[ERROR] CLEAN_infer_fasta.py not found at ${INFER_PY}" >&2
    exit 1
fi

if [ -d "${APP_DIR}/esm/.git" ]; then
    echo "[INFO] esm repo already cloned at ${APP_DIR}/esm"
else
    echo "[INFO] Cloning facebookresearch/esm into ${APP_DIR}/esm ..."
    git clone "${ESM_REPO_URL}" "${APP_DIR}/esm"
fi

(
    cd "${APP_DIR}/esm"
    CURRENT="$(git describe --tags --exact-match HEAD 2>/dev/null || echo unknown)"
    if [ "${CURRENT}" != "v1.0.2" ]; then
        echo "[INFO] Pinning esm to v1.0.2 (was ${CURRENT}) ..."
        git fetch --tags --quiet
        git checkout v1.0.2
    fi
)

cd "${APP_DIR}"
if [ -f "build.py" ]; then
    echo "[INFO] Building CLEAN package ..."
    python build.py install
fi

mkdir -p "${APP_DIR}/data/esm_data"
mkdir -p "${APP_DIR}/data/pretrained"
mkdir -p "${APP_DIR}/data/inputs"
mkdir -p "${APP_DIR}/results/inputs"

PRETRAINED_DIR="${APP_DIR}/data/pretrained"
SPLIT100_PTH="$(find "${PRETRAINED_DIR}" -maxdepth 2 -name 'split100*.pth' 2>/dev/null | head -1 || true)"

if [ -n "${SPLIT100_PTH}" ] && [ -s "${SPLIT100_PTH}" ]; then
    echo "[INFO] Pretrained weights already present: ${SPLIT100_PTH}"
else
    echo "[INFO] Downloading CLEAN pretrained weights from Google Drive ..."
    if ! command -v gdown &>/dev/null; then
        echo "[ERROR] gdown not found on PATH. Is the conda env active?" >&2
        exit 1
    fi
    cd "${PRETRAINED_DIR}"
    gdown "https://drive.google.com/uc?id=${WEIGHTS_GDRIVE_ID}" -O pretrained.zip
    if [ ! -s pretrained.zip ]; then
        echo "[ERROR] Pretrained weights download failed (empty file)." >&2
        echo "[ERROR] You may need to download manually from:" >&2
        echo "        https://drive.google.com/file/d/${WEIGHTS_GDRIVE_ID}/view" >&2
        echo "        and unzip its contents into ${PRETRAINED_DIR}" >&2
        rm -f pretrained.zip
        exit 1
    fi
    if command -v unzip &>/dev/null; then
        unzip -o pretrained.zip
    else
        python -c "import zipfile; zipfile.ZipFile('pretrained.zip').extractall()"
    fi
    rm -f pretrained.zip

    SUBDIR="$(find "${PRETRAINED_DIR}" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | head -1 || true)"
    if [ -n "${SUBDIR}" ] && [ -z "$(find "${PRETRAINED_DIR}" -maxdepth 1 -name '*.pth' 2>/dev/null)" ]; then
        echo "[INFO] Flattening sub-directory: ${SUBDIR}"
        mv "${SUBDIR}"/* "${PRETRAINED_DIR}"/ 2>/dev/null || true
        rmdir "${SUBDIR}" 2>/dev/null || true
    fi
fi

SPLIT100_PTH="$(find "${PRETRAINED_DIR}" -maxdepth 2 -name 'split100*.pth' 2>/dev/null | head -1 || true)"
if [ -z "${SPLIT100_PTH}" ] || [ ! -s "${SPLIT100_PTH}" ]; then
    echo "[ERROR] split100 weight file not found in ${PRETRAINED_DIR}" >&2
    echo "[ERROR] CLEAN inference will fail at runtime." >&2
    exit 1
fi
echo "[INFO] Found split100 weights: ${SPLIT100_PTH}"

touch "${INSTALL_DIR}/.setup_done"
echo "[INFO] CLEAN setup complete at ${INSTALL_DIR}"
