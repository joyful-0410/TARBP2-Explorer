import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test


# ==========================================
# TARBP2 Survival Analysis
# ==========================================

st.set_page_config(
    page_title="TARBP2 Survival Analysis",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 TARBP2 Survival Analysis")

st.markdown("""
本模組使用 TCGA Pan-Cancer RNA-seq 與臨床資料，
分析 TARBP2 表現量與癌症患者整體存活率（Overall Survival, OS）的關係。
""")

st.divider()


# ==========================================
# 1. Load and merge data
# ==========================================

@st.cache_data
def load_survival_data():
    expression_path = "data/TCGA/TARBP2_primary_tumor.csv"
    clinical_path = (
        "data/clinical/"
        "Survival_SupplementalTable_S1_20171025_xena_sp"
    )

    expression = pd.read_csv(expression_path)
    clinical = pd.read_csv(clinical_path, sep="\t")

    expression_required = [
        "sample",
        "TARBP2_expression",
        "_primary_disease"
    ]

    clinical_required = [
        "sample",
        "_PATIENT",
        "cancer type abbreviation",
        "OS",
        "OS.time"
    ]

    for col in expression_required:
        if col not in expression.columns:
            raise ValueError(f"Expression 缺少欄位：{col}")

    for col in clinical_required:
        if col not in clinical.columns:
            raise ValueError(f"Clinical 缺少欄位：{col}")

    expression["sample"] = expression["sample"].astype(str).str.strip()
    clinical["sample"] = clinical["sample"].astype(str).str.strip()
    clinical["_PATIENT"] = clinical["_PATIENT"].astype(str).str.strip()

    expression["TARBP2_expression"] = pd.to_numeric(
        expression["TARBP2_expression"], errors="coerce"
    )
    clinical["OS"] = pd.to_numeric(clinical["OS"], errors="coerce")
    clinical["OS.time"] = pd.to_numeric(clinical["OS.time"], errors="coerce")

    expression = expression.dropna(
        subset=["sample", "TARBP2_expression", "_primary_disease"]
    )
    clinical = clinical.dropna(
        subset=[
            "sample",
            "_PATIENT",
            "OS",
            "OS.time",
            "cancer type abbreviation"
        ]
    )

    clinical = clinical[
        clinical["OS"].isin([0, 1]) & (clinical["OS.time"] > 0)
    ].copy()

    expression = expression.drop_duplicates(subset=["sample"], keep="first")
    clinical = clinical.drop_duplicates(subset=["sample"], keep="first")

    merged = expression.merge(
        clinical,
        on="sample",
        how="inner",
        validate="one_to_one"
    )

    # 同一病患若有多筆樣本，只保留一筆；本資料使用 primary tumor expression。
    merged = merged.sort_values("sample")
    merged = merged.drop_duplicates(subset=["_PATIENT"], keep="first")

    if merged.empty:
        raise ValueError("Expression 與 Clinical 沒有成功配對的樣本。")

    return merged


try:
    df = load_survival_data()
except Exception as e:
    st.error(f"資料讀取或合併失敗：{e}")
    st.stop()


# ==========================================
# 2. Dataset overview
# ==========================================

st.subheader("1. Dataset Overview")

col1, col2, col3 = st.columns(3)

col1.metric("Matched Patients", f"{len(df):,}")
col2.metric("Death Events", f"{int(df['OS'].sum()):,}")
col3.metric("Cancer Types", int(df["cancer type abbreviation"].nunique()))

st.caption("分析單位為病患；同一病患若有多筆樣本，只保留一筆。")
st.divider()


# ==========================================
# 3. Select cancer
# ==========================================

st.subheader("2. Select Cancer Type")

cancer_options = sorted(
    df["cancer type abbreviation"].dropna().unique().tolist()
)

if not cancer_options:
    st.error("資料中沒有可供選擇的癌症種類。")
    st.stop()

# 若資料含 LUAD，預設選取 LUAD；否則預設第一個癌症種類。
default_index = cancer_options.index("LUAD") if "LUAD" in cancer_options else 0

selected_cancer = st.selectbox(
    "Select Cancer",
    cancer_options,
    index=default_index
)

cancer_df = df[
    df["cancer type abbreviation"] == selected_cancer
].copy()

if cancer_df.empty:
    st.warning("此癌症種類沒有可用資料。")
    st.stop()


# ==========================================
# 4. Prepare survival groups
# ==========================================

median_expression = cancer_df["TARBP2_expression"].median()

cancer_df["Expression Group"] = np.where(
    cancer_df["TARBP2_expression"] > median_expression,
    "High Expression",
    "Low Expression"
)

# 1 = High expression, 0 = Low expression
cancer_df["High_Group"] = (
    cancer_df["Expression Group"] == "High Expression"
).astype(int)

high_group = cancer_df[
    cancer_df["Expression Group"] == "High Expression"
].copy()

low_group = cancer_df[
    cancer_df["Expression Group"] == "Low Expression"
].copy()


# ==========================================
# 5. Group overview
# ==========================================

st.subheader("3. TARBP2 Expression Groups")

col1, col2, col3 = st.columns(3)

col1.metric("Total Patients", len(cancer_df))
col2.metric("High Expression", len(high_group))
col3.metric("Low Expression", len(low_group))

st.write(f"TARBP2 Median Expression: **{median_expression:.3f}**")

if high_group.empty or low_group.empty:
    st.error(
        "無法將此癌症分成高、低表現兩組（其中一組沒有樣本）。"
        "請檢查該癌症的 TARBP2 表現量分布。"
    )
    st.stop()

st.divider()


# ==========================================
# 6. Survival overview
# ==========================================

st.subheader("4. Survival Information")

summary = pd.DataFrame({
    "Group": ["High Expression", "Low Expression"],
    "Patients": [len(high_group), len(low_group)],
    "Death Events": [
        int(high_group["OS"].sum()),
        int(low_group["OS"].sum())
    ],
    "Median Follow-up (Days)": [
        high_group["OS.time"].median(),
        low_group["OS.time"].median()
    ]
})

st.dataframe(
    summary,
    use_container_width=True,
    hide_index=True
)


# ==========================================
# 7. Kaplan-Meier Survival Curve
# ==========================================

st.divider()
st.subheader("5. Kaplan-Meier Survival Analysis")

kmf_high = KaplanMeierFitter()
kmf_low = KaplanMeierFitter()

fig, ax = plt.subplots(figsize=(9, 6))

kmf_high.fit(
    durations=high_group["OS.time"],
    event_observed=high_group["OS"],
    label="High TARBP2 Expression"
)

kmf_low.fit(
    durations=low_group["OS.time"],
    event_observed=low_group["OS"],
    label="Low TARBP2 Expression"
)

kmf_high.plot_survival_function(ax=ax, ci_show=True)
kmf_low.plot_survival_function(ax=ax, ci_show=True)

ax.set_title(f"Kaplan-Meier Survival Curve - {selected_cancer}")
ax.set_xlabel("Time (Days)")
ax.set_ylabel("Overall Survival Probability")
ax.set_ylim(0, 1.05)
ax.grid(alpha=0.2)

plt.tight_layout()
st.pyplot(fig)
plt.close(fig)


# ==========================================
# 8. Log-rank test
# ==========================================

st.divider()
st.subheader("6. Log-rank Test")

logrank_result = logrank_test(
    high_group["OS.time"],
    low_group["OS.time"],
    event_observed_A=high_group["OS"],
    event_observed_B=low_group["OS"]
)

p_value = float(logrank_result.p_value)

col1, col2 = st.columns(2)

col1.metric("Log-rank P-value", f"{p_value:.4e}")

if p_value < 0.05:
    col2.success("Statistically Significant (P < 0.05)")
else:
    col2.info("Not Statistically Significant (P ≥ 0.05)")


# ==========================================
# 9. Cox Proportional Hazards Model
# ==========================================

st.divider()
st.subheader("7. Cox Proportional Hazards Regression")

cox_data = cancer_df[
    ["OS.time", "OS", "High_Group"]
].copy()

try:
    cph = CoxPHFitter()
    cph.fit(
        cox_data,
        duration_col="OS.time",
        event_col="OS"
    )

    cox_summary = cph.summary

    hr = float(cox_summary.loc["High_Group", "exp(coef)"])
    hr_lower = float(cox_summary.loc["High_Group", "exp(coef) lower 95%"])
    hr_upper = float(cox_summary.loc["High_Group", "exp(coef) upper 95%"])
    cox_p = float(cox_summary.loc["High_Group", "p"])

    col1, col2, col3 = st.columns(3)

    col1.metric("Hazard Ratio", f"{hr:.3f}")
    col2.metric("95% CI", f"{hr_lower:.3f} - {hr_upper:.3f}")
    col3.metric("Cox P-value", f"{cox_p:.4e}")

    st.dataframe(
        cox_summary,
        use_container_width=True
    )

    # 只有達統計顯著時才描述為顯著關聯。
    if cox_p < 0.05:
        if hr > 1:
            st.warning(
                "High TARBP2 expression is associated with a "
                "significantly higher estimated hazard of death."
            )
        else:
            st.success(
                "High TARBP2 expression is associated with a "
                "significantly lower estimated hazard of death."
            )
    else:
        if hr > 1:
            direction_text = (
                "The estimated hazard is higher in the high-expression group, "
                "but the difference is not statistically significant."
            )
        elif hr < 1:
            direction_text = (
                "The estimated hazard is lower in the high-expression group, "
                "but the difference is not statistically significant."
            )
        else:
            direction_text = "The estimated hazards are similar between groups."

        st.info(
            "No statistically significant association was detected "
            "between TARBP2 expression and overall survival. "
            + direction_text
        )

except Exception as e:
    st.warning(f"Cox regression could not be completed: {e}")


# ==========================================
# 10. Interpretation
# ==========================================

st.divider()
st.subheader("8. Interpretation")

if p_value < 0.05:
    st.success(
        "The Log-rank test detected a statistically significant "
        "survival difference between the two expression groups (P < 0.05)."
    )
else:
    st.info(
        "The Log-rank test did not detect a statistically significant "
        "survival difference between the two expression groups (P ≥ 0.05)."
    )

st.caption("""
Note:
1. High and low expression groups are divided by the median TARBP2 expression.
2. OS represents overall survival.
3. OS = 1 indicates death; OS = 0 indicates censored observations.
4. The Cox model in this version is univariate and is not adjusted for age,
   sex, or clinical stage.
5. Hazard ratio describes an association, not causation.
6. A non-significant result does not prove that no association exists.
""")


# ==========================================
# 11. Download
# ==========================================

st.divider()
st.subheader("9. Download Survival Data")

download_df = cancer_df.copy()

st.download_button(
    label="Download Survival Analysis Dataset",
    data=download_df.to_csv(index=False).encode("utf-8-sig"),
    file_name=f"TARBP2_{selected_cancer}_Survival.csv",
    mime="text/csv"
)

st.success("TARBP2 Survival Analysis completed.")
