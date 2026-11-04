#!/usr/bin/env bash

set -e

# 1. Asegurar que nos situamos en el directorio del proyecto raíz
cd "$(dirname "$0")"

# 2. Inicializar Conda para que permita usar "conda activate" dentro de un script de Bash
eval "$(conda shell.bash hook)"

# 3. Activar el entorno virtual local de Conda creado para la tarea
conda activate ./.venv

# 4. Ejecutar el script principal de Python de la carpeta src
python src/01_run_deseq2.py