# README.md - Pipeline de Análisis RNA-seq (PyDESeq2 + GSEApy)

- Estudiante: Ángela del Rocío De La Rosa Ruiz
- Lenguaje elegido: Python (versión 3.12)
- Repositorio GitHub: https://github.com/angderui/project_PyDESeq2_BAG_TCGA-BRCA.git

- Instrucciones:
1. Crear el entorno:
Abre una terminal en la raíz del proyecto y ejecuta:
``` bash
conda env create -p ./.venv -f environment.yml
```

2. Ejecutar el análisis:
Coloca los datos de entrada en data/raw/ (counts.csv.gz y metadata.csv.gz) y ejecuta el script principal main.sh en la terminal
``` bash
bash main.sh
```

3. Localizar los resultados:
 · Las tablas de resultados (TSV) y las figuras generadas (PNG) se encuentran en la carpeta results/.
 · Los registros de trazabilidad, logs y resúmenes de ejecución están en la carpeta logs/.

- Decisiones e incidencias relevantes:
    · Backend de Matplotlib: se configuró el backend "Agg" en el script de Python para permitir la generación y guardado automático de gráficos en formato PNG desde entornos de terminal sin interfaz gráfica.
    · Universo de fondo (Background): para el enriquecimiento funcional con GSEApy, se utilizó como universo de referencia el conjunto completo de genes analizados tras los filtros de calidad (res_df), garantizando así el rigor estadístico.
    · Entorno de ejecución: el flujo de trabajo se ha diseñado para ejecutarse de forma automatizada mediante un script de bash (main.sh), asegurando la reproducibilidad a partir del fichero "environment.yml".