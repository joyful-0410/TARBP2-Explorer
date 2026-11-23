import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import mannwhitneyu, wilcoxon


# ==========================================
# TARBP2 Tumor vs Normal Analysis
# ==========================================

st.set_page_config(
    page_title="TARBP2 Tumor vs Normal",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 TARBP2 Tumor vs Normal Analysis")

st.markdown("""
本模組利用 TCGA Pan-Cancer RNA-seq 資料，
比較不同癌症中 TARBP2 在腫瘤組織與正常組織的表現量差異。
""")

st.divider()


# ==========================================
# 1. Load data
# ==========================================

@st.cache_data
def load_data():

    file_path = "data/TCGA/TARBP2_TCGA_all_samples.csv"

    df = pd.read_csv(file_path)

    df["TARBP2_expression"] = pd.to_numeric(
        df["TARBP2_expression"],
        errors="coerce"
    )

    df = df.dropna(
        subset=[
            "TARBP2_expression",
            "_primary_disease",
            "sample_type"
        ]
    )

    return df


def add_patient_id(data):
    """
    TCGA sample barcode 的前 12 個字元代表病患層級 ID
    （例如 TCGA-AB-1234）。若沒有可辨識的 sample 欄位，
    則回傳 None，後續採獨立樣本比較並清楚標示限制。
    """
    sample_col = next(
        (col for col in ["sample", "sample_id", "Sample", "Sample ID"]
         if col in data.columns),
        None
    )
    if sample_col is None:
        # Keep the original dataframe; analysis will fall back to
        # independent-sample testing when no sample barcode is available.
        return data.copy()

    result = data.copy()
    result["Patient_ID"] = (
        result[sample_col].astype(str).str.strip().str[:12]
    )
    invalid = (
        result["Patient_ID"].isna()
        | result["Patient_ID"].isin(["", "nan", "None"])
    )
    result.loc[invalid, "Patient_ID"] = np.nan
    return result


try:

    df = load_data()
    df = add_patient_id(df)

except FileNotFoundError:

    st.error("找不到 TCGA 資料，請先執行 merge_tcga.py")

    st.stop()


# ==========================================
# 2. Select sample types
# ==========================================

tumor = df[
    df["sample_type"] == "Primary Tumor"
].copy()

normal = df[
    df["sample_type"] == "Solid Tissue Normal"
].copy()


# ==========================================
# 3. Overview
# ==========================================

st.subheader("1. Dataset Overview")

col1, col2, col3 = st.columns(3)

col1.metric(
    "Primary Tumor Samples",
    f"{len(tumor):,}"
)

col2.metric(
    "Normal Tissue Samples",
    f"{len(normal):,}"
)

col3.metric(
    "Gene",
    "TARBP2"
)

st.divider()


# ==========================================
# 4. Cancer-wise statistical analysis
# ==========================================

st.subheader("2. Cancer-wise Statistical Comparison")


def calculate_statistics(data):
    results = []

    cancers = sorted(data["_primary_disease"].dropna().unique())

    for cancer in cancers:
        cancer_data = data[data["_primary_disease"] == cancer].copy()

        tumor_data = cancer_data[
            cancer_data["sample_type"] == "Primary Tumor"
        ].copy()
        normal_data = cancer_data[
            cancer_data["sample_type"] == "Solid Tissue Normal"
        ].copy()

        tumor_n = len(tumor_data)
        normal_n = len(normal_data)

        tumor_values = tumor_data["TARBP2_expression"].dropna()
        normal_values = normal_data["TARBP2_expression"].dropna()

        if tumor_n < 5 or normal_n < 5:
            continue

        matched_pairs = 0
        test_method = "Mann–Whitney U（獨立樣本）"
        analysis_note = (
            "未使用配對檢定；結果比較兩組樣本分布，"
            "不代表病患內部變化。"
        )

        # If TCGA sample barcodes are available, identify patients with
        # both primary-tumor and solid-tissue-normal samples.
        paired_p_value = np.nan
        paired_tumor = pd.Series(dtype=float)
        paired_normal = pd.Series(dtype=float)

        if "Patient_ID" in cancer_data.columns:
            tumor_patient = (
                tumor_data.dropna(subset=["Patient_ID"])
                .groupby("Patient_ID")["TARBP2_expression"]
                .median()
            )
            normal_patient = (
                normal_data.dropna(subset=["Patient_ID"])
                .groupby("Patient_ID")["TARBP2_expression"]
                .median()
            )

            common_patients = tumor_patient.index.intersection(
                normal_patient.index
            )
            matched_pairs = len(common_patients)

            if matched_pairs >= 5:
                paired_tumor = tumor_patient.loc[common_patients]
                paired_normal = normal_patient.loc[common_patients]
                differences = paired_tumor - paired_normal

                # Wilcoxon signed-rank test is appropriate for paired
                # observations. Handle the all-zero-difference case.
                if np.allclose(differences.to_numpy(), 0):
                    statistic, p_value = 0.0, 1.0
                else:
                    statistic, p_value = wilcoxon(
                        paired_tumor,
                        paired_normal,
                        alternative="two-sided"
                    )

                tumor_median = paired_tumor.median()
                normal_median = paired_normal.median()
                median_difference = differences.median()
                test_method = "Wilcoxon signed-rank（配對樣本）"
                analysis_note = (
                    "以同一病患的腫瘤與正常樣本配對比較；"
                    "每位病患若有多筆同類樣本，先取表現量中位數。"
                )
                tumor_n_used = matched_pairs
                normal_n_used = matched_pairs
            else:
                statistic, p_value = mannwhitneyu(
                    tumor_values,
                    normal_values,
                    alternative="two-sided"
                )
                tumor_median = tumor_values.median()
                normal_median = normal_values.median()
                median_difference = tumor_median - normal_median
                tumor_n_used = len(tumor_values)
                normal_n_used = len(normal_values)
        else:
            statistic, p_value = mannwhitneyu(
                tumor_values,
                normal_values,
                alternative="two-sided"
            )
            tumor_median = tumor_values.median()
            normal_median = normal_values.median()
            median_difference = tumor_median - normal_median
            tumor_n_used = len(tumor_values)
            normal_n_used = len(normal_values)
            analysis_note = (
                "資料中找不到可辨識的 TCGA sample barcode，"
                "因此採獨立樣本比較；無法確認或利用病患配對關係。"
            )

        results.append({
            "Cancer Type": cancer,
            "Tumor N": tumor_n,
            "Normal N": normal_n,
            "Matched Pairs": matched_pairs,
            "Tumor N Used": tumor_n_used,
            "Normal N Used": normal_n_used,
            "Tumor Median": tumor_median,
            "Normal Median": normal_median,
            "Median Difference": median_difference,
            "Test Method": test_method,
            "P-value": p_value,
            "Analysis Note": analysis_note
        })

    return pd.DataFrame(results)


# Run the cancer-wise analysis before applying multiple-testing correction.
results = calculate_statistics(df)


# ==========================================
# 5. Benjamini-Hochberg correction
# ==========================================

def benjamini_hochberg(p_values):

    p_values = np.asarray(p_values, dtype=float)

    n = len(p_values)

    if n == 0:
        return np.array([])

    order = np.argsort(p_values)

    ranked_p = p_values[order]

    adjusted = ranked_p * n / np.arange(1, n + 1)

    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]

    adjusted = np.clip(adjusted, 0, 1)

    q_values = np.empty(n)

    q_values[order] = adjusted

    return q_values


if not results.empty:

    results["Adjusted P-value"] = benjamini_hochberg(
        results["P-value"]
    )

    results["Direction"] = np.where(
        results["Median Difference"] > 0,
        "Higher in Tumor",
        np.where(
            results["Median Difference"] < 0,
            "Lower in Tumor",
            "No Median Difference"
        )
    )

    results = results.sort_values(
        "Adjusted P-value"
    )


# ==========================================
# 6. Display statistical results
# ==========================================

if results.empty:

    st.warning(
        "目前沒有足夠的癌症樣本同時符合腫瘤與正常組織各至少 5 個樣本的條件。"
    )

    st.stop()


st.write(
    f"符合分析條件的癌症種類：{len(results)}"
)

st.caption(
    "統計方法：若同一病患同時有腫瘤與正常樣本且配對數至少 5，"
    "使用 Wilcoxon signed-rank test；否則使用 Mann–Whitney U test。"
    "所有癌症種類的 P-value 再以 Benjamini–Hochberg 方法進行 FDR 校正。"
)

st.dataframe(
    results,
    use_container_width=True,
    hide_index=True
)


# ==========================================
# 7. Significant results
# ==========================================

st.divider()

st.subheader("3. Statistical Significance")

significant = results[
    results["Adjusted P-value"] < 0.05
]

col1, col2 = st.columns(2)

col1.metric(
    "Analyzed Cancer Types",
    len(results)
)

col2.metric(
    "FDR < 0.05",
    len(significant)
)

if len(significant) > 0:

    st.success(
        f"共有 {len(significant)} 種癌症在多重比較校正後仍達到 FDR < 0.05。"
    )

else:

    st.info(
        "目前沒有癌症在多重比較校正後達到 FDR < 0.05。"
    )


# ==========================================
# 8. Select cancer
# ==========================================

st.divider()

st.subheader("4. Individual Cancer Comparison")

cancer_options = results["Cancer Type"].tolist()

selected_cancer = st.selectbox(
    "Select Cancer Type",
    cancer_options
)


selected_data = df[
    df["_primary_disease"] == selected_cancer
]


tumor_values = selected_data.loc[
    selected_data["sample_type"] == "Primary Tumor",
    "TARBP2_expression"
].dropna()


normal_values = selected_data.loc[
    selected_data["sample_type"] == "Solid Tissue Normal",
    "TARBP2_expression"
].dropna()


selected_result = results[
    results["Cancer Type"] == selected_cancer
].iloc[0]


# ==========================================
# 9. Metrics
# ==========================================

col1, col2, col3 = st.columns(3)

col1.metric(
    "Tumor Median",
    f"{tumor_values.median():.3f}"
)

col2.metric(
    "Normal Median",
    f"{normal_values.median():.3f}"
)

col3.metric(
    "Adjusted P-value",
    f"{selected_result['Adjusted P-value']:.4g}"
)


st.write(
    f"Tumor samples: {len(tumor_values)}"
)

st.write(
    f"Normal samples: {len(normal_values)}"
)


# ==========================================
# 10. Box Plot
# ==========================================

fig, ax = plt.subplots(figsize=(8, 5))

ax.boxplot(
    [
        tumor_values,
        normal_values
    ],
    tick_labels=[
        "Primary Tumor",
        "Solid Tissue Normal"
    ],
    showfliers=False
)

ax.set_title(
    f"TARBP2 Expression in {selected_cancer}"
)

ax.set_ylabel("TARBP2 Expression")

ax.set_xlabel("Sample Type")

plt.tight_layout()

st.pyplot(fig)

plt.close(fig)


# ==========================================
# 11. Interpretation
# ==========================================

st.divider()

st.subheader("5. Statistical Interpretation")

difference = selected_result["Median Difference"]

if difference > 0:

    st.write(
        "TARBP2 median expression is higher in primary tumor samples."
    )

elif difference < 0:

    st.write(
        "TARBP2 median expression is lower in primary tumor samples."
    )

else:

    st.write(
        "Tumor and normal samples have the same median expression."
    )


if selected_result["Adjusted P-value"] < 0.05:

    st.success(
        "The difference remains statistically significant "
        "after Benjamini-Hochberg correction (FDR < 0.05)."
    )

else:

    st.info(
        "The difference does not reach FDR < 0.05 "
        "after multiple-testing correction."
    )


st.caption(
    "Note: This is an exploratory analysis. When sufficient matched "
    "TCGA patient pairs are identifiable, paired samples are analyzed "
    "with the Wilcoxon signed-rank test; otherwise, groups are compared "
    "as independent samples with the Mann–Whitney U test. "
    "Clinical covariates are not adjusted, and statistical significance "
    "does not establish biological causality."
)


# ==========================================
# 12. Download
# ==========================================

st.divider()

st.subheader("6. Download Results")

st.download_button(
    label="Download Tumor vs Normal Statistics",

    data=results.to_csv(
        index=False
    ).encode("utf-8-sig"),

    file_name="TARBP2_Tumor_Normal_Statistics.csv",

    mime="text/csv"
)

st.success(
    "TARBP2 Tumor vs Normal Analysis completed."
)