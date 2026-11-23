import pandas as pd
import os

# =====================================
# 1. 讀取資料
# =====================================

expression_file = "data/expression/TARBP2_candidates.csv"

phenotype_file = (
    "data/expression/"
    "TCGA_phenotype_denseDataOnlyDownload.tsv"
)

print("正在讀取 TARBP2 表現量...")
expression = pd.read_csv(expression_file)

print("正在讀取 TCGA Phenotype...")
phenotype = pd.read_csv(phenotype_file, sep="\t")


# =====================================
# 2. 整理 TARBP2 表現量
# =====================================

# 第一欄是基因名稱，不是樣本名稱
expression = expression.rename(
    columns={expression.columns[0]: "gene_id"}
)

gene_name = expression.loc[0, "gene_id"]

# 將寬格式轉換成直式格式
expression_long = expression.melt(
    id_vars="gene_id",
    var_name="sample",
    value_name="TARBP2_expression"
)

expression_long["TARBP2_expression"] = pd.to_numeric(
    expression_long["TARBP2_expression"],
    errors="coerce"
)

print("基因名稱：", gene_name)
print("表現量樣本數：", len(expression_long))


# =====================================
# 3. 整理 Phenotype
# =====================================

phenotype = phenotype[
    [
        "sample",
        "sample_type_id",
        "sample_type",
        "_primary_disease"
    ]
].copy()

phenotype["sample"] = phenotype["sample"].astype(str).str.strip()

# 移除重複樣本
phenotype = phenotype.drop_duplicates(
    subset="sample"
)


# =====================================
# 4. 合併資料
# =====================================

expression_long["sample"] = (
    expression_long["sample"].astype(str).str.strip()
)

merged = expression_long.merge(
    phenotype,
    on="sample",
    how="inner"
)

print("\n成功配對數：", len(merged))


# =====================================
# 5. 排除缺失資料
# =====================================

merged = merged.dropna(
    subset=[
        "TARBP2_expression",
        "_primary_disease",
        "sample_type"
    ]
)


# =====================================
# 6. 分離樣本類型
# =====================================

primary_tumor = merged[
    merged["sample_type"] == "Primary Tumor"
].copy()

normal = merged[
    merged["sample_type"] == "Solid Tissue Normal"
].copy()

metastatic = merged[
    merged["sample_type"] == "Metastatic"
].copy()


# =====================================
# 7. 建立輸出資料夾
# =====================================

output_dir = "data/TCGA"
os.makedirs(output_dir, exist_ok=True)


# =====================================
# 8. 輸出資料
# =====================================

merged.to_csv(
    f"{output_dir}/TARBP2_TCGA_all_samples.csv",
    index=False
)

primary_tumor.to_csv(
    f"{output_dir}/TARBP2_primary_tumor.csv",
    index=False
)

normal.to_csv(
    f"{output_dir}/TARBP2_normal_tissue.csv",
    index=False
)

metastatic.to_csv(
    f"{output_dir}/TARBP2_metastatic.csv",
    index=False
)


# =====================================
# 9. 癌症表現量統計
# =====================================

summary = (
    primary_tumor
    .groupby("_primary_disease")
    ["TARBP2_expression"]
    .agg(
        sample_count="count",
        mean="mean",
        median="median",
        std="std"
    )
    .sort_values("median", ascending=False)
)

summary.to_csv(
    f"{output_dir}/TARBP2_cancer_summary.csv"
)


# =====================================
# 10. 顯示結果
# =====================================

print("\n========== 分析完成 ==========")

print("全部配對樣本：", len(merged))
print("原發性腫瘤：", len(primary_tumor))
print("正常組織：", len(normal))
print("轉移性腫瘤：", len(metastatic))

print(
    "\n原發性腫瘤癌症種類：",
    primary_tumor["_primary_disease"].nunique()
)

print("\n各癌症表現量統計：")
print(summary.to_string())

print("\n輸出檔案位置：data/TCGA/")