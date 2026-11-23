import pandas as pd

# 讀取 TARBP2 表現量資料
expression = pd.read_csv(
    "data/expression/TARBP2_candidates.csv"
)

# 讀取 TCGA Phenotype
phenotype = pd.read_csv(
    "data/expression/TCGA_phenotype_denseDataOnlyDownload.tsv",
    sep="\t"
)

# 取得表現量資料中的樣本名稱
expression_samples = set(expression.columns[1:])

# 取得 Phenotype 樣本名稱
phenotype_samples = set(phenotype["sample"].dropna())

# 計算交集
matched = expression_samples & phenotype_samples

print("=" * 45)
print("TCGA 樣本配對檢查")
print("=" * 45)

print("TARBP2 樣本數：", len(expression_samples))
print("Phenotype 樣本數：", len(phenotype_samples))
print("成功配對數：", len(matched))

# 查看 Phenotype 樣本類型
print("\n樣本類型分布：")

print(
    phenotype["sample_type"]
    .value_counts()
    .to_string()
)

# 查看疾病種類數量
print("\n癌症種類數：")

print(
    phenotype["_primary_disease"].nunique()
)

# 查看成功配對的樣本類型
matched_data = phenotype[
    phenotype["sample"].isin(matched)
]

print("\n成功配對的樣本類型：")

print(
    matched_data["sample_type"]
    .value_counts()
    .to_string()
)

print("\n前10筆成功配對資料：")

print(
    matched_data.head(10).to_string(index=False)
)