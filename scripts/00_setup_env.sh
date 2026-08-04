#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

if ! command -v micromamba >/dev/null 2>&1; then
  mkdir -p .tools/bin .tools/micromamba
  if [ ! -x .tools/bin/micromamba ]; then
    echo "Installing micromamba locally under .tools/bin"
    if command -v curl >/dev/null 2>&1; then
      curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C .tools/micromamba bin/micromamba
    elif command -v wget >/dev/null 2>&1; then
      wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C .tools/micromamba bin/micromamba
    else
      echo "Neither curl nor wget was found. Install one of them first, then rerun this script."
      exit 1
    fi
    cp .tools/micromamba/bin/micromamba .tools/bin/micromamba
    chmod +x .tools/bin/micromamba
  fi
  export PATH="$PWD/.tools/bin:$PATH"
fi

micromamba create -y -p .mamba_vina -f environment-vina.yml

echo "Environment ready. Use:"
echo "micromamba run -p .mamba_vina python scripts/01_check_inputs.py"
