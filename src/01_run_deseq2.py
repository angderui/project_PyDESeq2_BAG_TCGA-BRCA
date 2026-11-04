#!/usr/bin/env python3

"""
Script de análisis diferencial de expresión con PyDESeq2 y enriquecimiento funcional con GSEApy de tumor frente a tejido 
normal con los datos TCGA-BRCA
Asignatura: Bioinformática y Análisis Genómico - Grado en Bioquímica
"""
import sys
from pathlib import Path
import importlib.metadata
import numpy as np
import pandas as pd

# Configurar matplotlib en modo no interactivo (Headless/Agg) para guardar figuras sin interfaz gráfica
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Importar librerías de análisis diferencial y enriquecimiento
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp

def main():
    # 0. PREPARACIÓN Y PARÁMETROS
    print("[0/8] Preparando directorios y parámetros")
    data_dir = Path("data/raw")
    results_dir = Path("results")
    logs_dir = Path("logs")

    results_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    counts_path = data_dir / "counts.csv.gz"
    metadata_path = data_dir / "metadata.csv.gz"

    condition_column = "site"
    plate_column = "plate"
    condition_test = "Primary_Solid_Tumor"
    condition_ref = "Solid_Tissue_Normal"

    min_count = 10
    min_samples = 20
    alpha = 0.05
    lfc_threshold = 1.0

    # 1. LEER DATOS DE ENTRADA
    print("[1/8] Leyendo matrices de conteos y metadatos (TCGA-BRCA)")
    counts_df = pd.read_csv(counts_path, index_col=0, low_memory=False)
    metadata = pd.read_csv(metadata_path, index_col=0, low_memory=False)

    print(f"Dimensiones iniciales de counts: {counts_df.shape} (Esperado: ~20502 genes x 75 muestras)")

    # 2. ALINEAR MUESTRAS
    print("[2/8] Alineando muestras entre conteos y metadatos")

    # Reordenar metadata para que coincida exactamente con las columnas de counts
    common_samples = counts_df.columns.intersection(metadata.index)
    counts_df = counts_df.loc[:, common_samples]
    metadata = metadata.loc[common_samples]

    # Convertir variables a categóricas
    metadata[condition_column] = metadata[condition_column].astype("category")
    metadata[plate_column] = metadata[plate_column].astype("category")

    # Asegurar nivel de referencia
    metadata[condition_column] = metadata[condition_column].cat.set_categories([condition_ref, condition_test], 
    ordered=True)

    assert list(counts_df.columns) == list(metadata.index), "Error: Las muestras no coinciden."

    # 3. FILTRADO Y OBJETO DESeq2
    print("[3/8] Filtrando genes de baja cuenta y creando objeto PyDESeq2")
    
    # PyDESeq2 requiere que las muestras estén en las filas y los genes en las columnas
    counts_mat = counts_df.T

    # Filtrado: gen con al menos min_count lecturas en al menos min_samples muestras
    keep_genes = (counts_df >= min_count).sum(axis=1) >= min_samples
    counts_filtered_df = counts_df.loc[keep_genes]
    counts_mat_filtered = counts_filtered_df.T

    print(f"Genes tras el filtrado: {counts_filtered_df.shape[0]} (Esperado: 17606)")

    # Crear el objeto DeseqDataSet indicando el diseño experimental
    dds = DeseqDataSet(counts=counts_mat_filtered, metadata=metadata, design_factors=[condition_column, plate_column], 
    refit_cooks=True, n_cpus=4)

    # 4. EJECUTAR DESeq
    print("[4/8] Ejecutando el modelo DESeq2")
    dds.deseq2()

    # Guardar matriz de conteos normalizados
    norm_counts = pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names)
    norm_counts.to_csv(results_dir / "normalized_counts.tsv", sep="\t")

    # Generar y guardar el gráfico de dispersión
    plt.figure(figsize=(8, 6))
    dds.plot_dispersions()
    plt.title("Dispersion Estimates")
    plt.tight_layout()
    plt.savefig(results_dir / "plot_dispersion.png", dpi=300)
    plt.close()

    # 5. RESULTADOS Y TABLAS
    print("[5/8] Extrayendo resultados, aplicando contrastes y generando tablas")

    # Configurar contraste estadístico: Tumor vs Normal
    stat_res = DeseqStats(dds, contrast = [condition_column, condition_test, condition_ref],
    alpha=alpha, n_cpus=4)
    stat_res.summary()

    res_df = stat_res.results_df.copy()
    res_df["gene_id"] = res_df.index

    # Cálculo de columna "minus_log10_padj" para futuras gráficas
    res_df["minus_log10_padj"] = -np.log10(res_df["padj"] + 1e-300)

    # Guardar resultados completos
    res_df.to_csv(results_dir / "deseq2_all_results.tsv", sep="\t", index=False)

    # Ordenar por padj (los valores NaN se colocan al final)
    res_sorted = res_df.sort_values(by="padj", na_position="last")
    res_sorted.to_csv(results_dir / "deseq2_sig_padj_0.05.tsv", sep="\t", index=False)

    # Filtrar genes significativos (sig_df) y guardarlos
    sig_df = res_sorted[res_sorted["padj"] < alpha].copy()
    sig_df.to_csv(results_dir / "deseq2_sig_genes.tsv", sep="\t", index=False)

    # Filtrar genes diferencialmente expresados (de_df) y guardarlos
    de_df = sig_df[sig_df["log2FoldChange"].abs() > lfc_threshold].copy()

    # Separar genes Up y Down regulados
    up_df = de_df[de_df["log2FoldChange"] > 0].sort_values(by="padj")
    down_df = de_df[de_df["log2FoldChange"] < 0].sort_values(by="padj")

    # Guardar tablas Up y Down completas en TSV
    up_df.to_csv(results_dir / "deseq2_up_genes.tsv", sep="\t", index=False)
    down_df.to_csv(results_dir / "deseq2_down_genes.tsv", sep="\t", index=False)

    # Extraer los 10 primeros genes de cada dirección
    top10_up = up_df.head(10).copy()
    top10_down = down_df.head(10).copy()

    top10_up.to_csv(results_dir / "top10_up_genes.tsv", sep="\t", index=False)
    top10_down.to_csv(results_dir / "top10_down_genes.tsv", sep="\t", index=False)

    # Añadir la columna 'top10_label' a de_df
    de_df["top10_label"] = "no"
    de_df.loc[top10_up.index, "top10_label"] = "up"
    de_df.loc[top10_down.index, "top10_label"] = "down"

    # Guardar de_df actualizado con su nueva columna en TSV
    de_df.to_csv(results_dir / "deseq2_de_padj_0.05_1fc_1.tsv", sep="\t", index=False)

    # 6. GRÁFICOS
    print("[6/8] Generando gráficos (Top 10, Volcano y MA plot)")

    # Gráficos de barras horizontales para el Top 10 de genes Up y Down
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    top10_up_sorted = top10_up.sort_values(by="log2FoldChange", ascending=True)
    axes[0].barh(top10_up_sorted["gene_id"], top10_up_sorted["log2FoldChange"], color="firebrick")
    axes[0].set_title("Top 10 Up (Tumor > Normal)")
    axes[0].set_xlabel("log2FoldChange")

    top10_down_sorted = top10_down.sort_values(by="log2FoldChange", ascending=False)
    axes[1].barh(top10_down_sorted["gene_id"], top10_down_sorted["log2FoldChange"], color="steelblue")
    axes[1].set_title("Top 10 Down (Tumor > Normal)")
    axes[1].set_xlabel("log2FoldChange")

    plt.tight_layout()
    plt.savefig(results_dir / "plot_top10_genes.png", dpi=300)
    plt.close()

    # Preparar datos válidos para Volcano y MA plots
    res_plot = res_df.dropna(subset=["padj", "log2FoldChange"]).copy()

    colors = np.where((res_plot["padj"] < alpha) & (res_plot["log2FoldChange"] > lfc_threshold), "firebrick", 
    np.where((res_plot["padj"] < alpha) & (res_plot["log2FoldChange"] < -lfc_threshold), "steelblue", "darkgray"))

    # Volcano Plot
    plt.figure(figsize=(8, 6))

    # 1. Definir máscaras para separar los grupos según significación y fold change
    up_mask = (res_plot["padj"] < alpha) & (res_plot["log2FoldChange"] > lfc_threshold)
    down_mask = (res_plot["padj"] < alpha) & (res_plot["log2FoldChange"] < -lfc_threshold)
    ns_mask = ~(up_mask | down_mask)  # No significativo

    # 2. Dibujar cada grupo por separado con su respectivo 'label' para la leyenda
    plt.scatter(res_plot.loc[ns_mask, "log2FoldChange"], res_plot.loc[ns_mask, "minus_log10_padj"], 
    c="darkgray", s=10, alpha=0.5, label="No significativo")
    plt.scatter(res_plot.loc[up_mask, "log2FoldChange"], res_plot.loc[up_mask, "minus_log10_padj"], 
    c="firebrick", s=10, alpha=0.7, label="Up")
    plt.scatter(res_plot.loc[down_mask, "log2FoldChange"], res_plot.loc[down_mask, "minus_log10_padj"], 
    c="steelblue", s=10, alpha=0.7, label="Down")

    # 3. Líneas de umbral
    plt.axvline(lfc_threshold, color="black", linestyle="--", alpha=0.05)
    plt.axvline(-lfc_threshold, color="black", linestyle="--", alpha=0.05)
    plt.axhline(-np.log10(alpha), color="black", linestyle="--", alpha=0.05)

    # 4. Etiquetas de texto para los Top 10 Up y Down
    top10_volcano = pd.concat([top10_up, top10_down])

    for idx, row in top10_volcano.iterrows():
        plt.text(x=row["log2FoldChange"], y=row["minus_log10_padj"], s=row["gene_id"], fontsize=8, weight='bold', 
        ha='right', va='bottom')

    plt.xlabel("log2FoldChange")
    plt.ylabel("-log10(padj)")
    plt.title("Volcano Plot (Tumor vs Normal)")
    
    # 5. Añadir la leyenda al gráfico
    plt.legend(loc="upper left", frameon=True)

    plt.tight_layout()
    plt.savefig(results_dir / "plot_volcano.png", dpi=300)
    plt.close()

    # MA Plot
    plt.figure(figsize=(8, 6))
    plt.scatter(res_plot["baseMean"], res_plot["log2FoldChange"], c=colors, s=10, alpha=0.7)
    plt.axhline(0, color="black", linewidth=1)
    plt.xscale("log")
    plt.xlabel("mean of normalized counts")
    plt.ylabel("log2 fold change")
    plt.title("DESeq2 MA Plot")
    plt.tight_layout()
    plt.savefig(results_dir / "plot_ma.png", dpi=300)
    plt.close()

    # 7. GEN DESTACADO
    print("[7/8] Guardando información del gen destacado")

    top_gene = de_df.iloc[0]
    with open(results_dir / "top_gene_context.txt", "w") as f:
        f.write(f"Gen destacado: {top_gene['gene_id']}\n")
        f.write(f"log2FoldChange: {top_gene['log2FoldChange']}\n")
        f.write(f"padj: {top_gene['padj']}\n")

    # 8. ENRIQUECIMIENTO GO Y REACTOME
    print("[8/8] Ejecutando enriquecimiento funcional (GO y Reactome)")

    def run_gsea_enrichment(gene_list, library_name, out_prefix):
        if len(gene_list) == 0:
            return
        try:
            # Consultar Enrichr usando el universo completo de genes (background)
            enr = gp.enrichr(gene_list=list(gene_list), gene_sets=library_name, organism="human",
            background=list(res_df["gene_id"]), outdir=None, no_plot=True)
            df_res = enr.res2d

            # Exportar la tabla de resultados completa a formato TSV
            df_res.to_csv(results_dir / f"{out_prefix}.tsv", sep="\t", index=False)

            if not df_res.empty:
                # Generar y guardar el gráfico dotplot con el Top 10 de términos
                gp.dotplot(df_res, title=out_prefix, cmap="viridis", max_term=10, 
                ofname=(results_dir / f"plot_{out_prefix}.png"))
                plt.close()
        except Exception as e:
            print(f"Aviso en enriquecimiento {out_prefix}: {e}")

    # Preparar las listas de genes diferencialmente expresados sin valores NaN
    up_genes_list = up_df["gene_id"].dropna()
    down_genes_list = down_df["gene_id"].dropna()

    # Ejecutar la función para las 4 combinaciones requeridas (GO, Reactome, Up y Down)
    run_gsea_enrichment(up_genes_list, "GO_Biological_Process_2023", "go_bp_up")
    run_gsea_enrichment(down_genes_list, "GO_Biological_Process_2023", "go_bp_down")
    run_gsea_enrichment(up_genes_list, "Reactome_2022", "reactome_up")
    run_gsea_enrichment(down_genes_list, "Reactome_2022", "reactome_down")

    # 9. TRAZABILIDAD
    print("Guardando registros de trazabilidad y logs")

    # Registrar las versiones de Python y librerías
    versions_info = f"""Python version: {sys.version}
    pandas: {importlib.metadata.version('pandas')}
    numpy: {importlib.metadata.version('numpy')}
    pydeseq2: {importlib.metadata.version('pydeseq2')}
    gseapy: {importlib.metadata.version('gseapy')}
    matplotlib: {importlib.metadata.version('matplotlib')}"""

    with open(logs_dir / "sessionInfo.txt", "w") as f:
        f.write(versions_info)
    
    # Guardar un resumen estadístico con los recuentos clave de cada etapa del pipeline
    summary_info = f"""=== RUN SUMMARY ===
    Total genes input: {counts_df.shape[0]}
    Filtered genes: {counts_filtered_df.shape[0]}
    Significant genes (padj < 0.05): {len(sig_df)}
    Differentially Expressed (padj < 0.05, |log2FC| > 1): {len(de_df)}
    Upregulated: {len(up_df)}
    Downregulated: {len(down_df)}"""

    with open(logs_dir / "run_summary.txt", "w") as f:
        f.write(summary_info)

    print("Pipeline ejecutado con éxito. Todos los resultados se encuentran en results/ y logs/.")


if __name__ == "__main__":
    main()