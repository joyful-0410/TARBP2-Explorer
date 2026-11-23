import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test


# =========================================================
# TARBP2 KIRC Deep-Dive Survival Analysis
# =========================================================

st.set_page_config(
    page_title="TARBP2 KIRC Deep-Dive",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 TARBP2 in KIRC: Survival Analysis")
st.markdown("""
本模組聚焦於 KIRC（Kidney Renal Clear Cell Carcinoma，腎透明細胞癌），
進一步檢視 TARBP2 表現量與 Overall Survival（OS）的關聯。

分析包含 Kaplan-Meier 曲線、Log-rank test、單變項 Cox、
臨床變項校正後的多變項 Cox，以及比例風險假設檢查。
""")
st.info(
    "注意：這是 TCGA 觀察性資料分析。統計關聯不代表因果關係；"
    "多變項分析結果也會受到臨床資料缺失及變項定義影響。"
)

st.divider()


# =========================================================
# 1. Load and merge data
# =========================================================

@st.cache_data
def load_data():
    expression_path = "data/TCGA/TARBP2_primary_tumor.csv"
    clinical_path = (
        "data/clinical/"
        "Survival_SupplementalTable_S1_20171025_xena_sp"
    )

    expression = pd.read_csv(expression_path)
    clinical = pd.read_csv(clinical_path, sep="\t")

    expression_required = [
        "sample", "TARBP2_expression", "_primary_disease"
    ]
    clinical_required = [
        "sample", "_PATIENT", "cancer type abbreviation",
        "OS", "OS.time"
    ]

    for col in expression_required:
        if col not in expression.columns:
            raise ValueError(f"Expression 資料缺少欄位：{col}")
    for col in clinical_required:
        if col not in clinical.columns:
            raise ValueError(f"Clinical 資料缺少欄位：{col}")

    expression["sample"] = expression["sample"].astype(str).str.strip()
    clinical["sample"] = clinical["sample"].astype(str).str.strip()
    clinical["_PATIENT"] = clinical["_PATIENT"].astype(str).str.strip()

    expression["TARBP2_expression"] = pd.to_numeric(
        expression["TARBP2_expression"], errors="coerce"
    )
    clinical["OS"] = pd.to_numeric(clinical["OS"], errors="coerce")
    clinical["OS.time"] = pd.to_numeric(
        clinical["OS.time"], errors="coerce"
    )

    expression = expression.dropna(
        subset=["sample", "TARBP2_expression"]
    ).drop_duplicates(subset=["sample"], keep="first")

    clinical = clinical.dropna(
        subset=[
            "sample", "_PATIENT", "cancer type abbreviation",
            "OS", "OS.time"
        ]
    )
    clinical = clinical[
        clinical["OS"].isin([0, 1]) & (clinical["OS.time"] > 0)
    ].drop_duplicates(subset=["sample"], keep="first")

    merged = expression.merge(
        clinical, on="sample", how="inner", validate="one_to_one"
    )

    # Restrict to KIRC.
    kirc = merged[
        merged["cancer type abbreviation"].astype(str).str.upper() == "KIRC"
    ].copy()

    if kirc.empty:
        raise ValueError(
            "沒有找到 KIRC 配對資料。請確認 cancer type abbreviation 欄位及資料檔。"
        )

    # The expression file is already restricted to primary-tumor samples.
    # If a patient has multiple matched primary-tumor samples, summarize
    # TARBP2 expression by the patient-level median instead of arbitrarily
    # keeping the alphabetically first sample.
    aggregation = {}
    for col in kirc.columns:
        if col == "_PATIENT":
            continue
        elif col == "TARBP2_expression":
            aggregation[col] = "median"
        else:
            aggregation[col] = "first"

    kirc = (
        kirc.sort_values("sample")
        .groupby("_PATIENT", as_index=False, dropna=False)
        .agg(aggregation)
    )

    if kirc["_PATIENT"].duplicated().any():
        raise ValueError("病患層級整理後仍存在重複 _PATIENT，請檢查資料。")

    return kirc


try:
    df = load_data()
except Exception as e:
    st.error(f"資料讀取或合併失敗：{e}")
    st.stop()


# =========================================================
# 2. Prepare expression groups and clinical covariates
# =========================================================

df["TARBP2_expression"] = pd.to_numeric(
    df["TARBP2_expression"], errors="coerce"
)
df["OS"] = pd.to_numeric(df["OS"], errors="coerce")
df["OS.time"] = pd.to_numeric(df["OS.time"], errors="coerce")

df = df.dropna(
    subset=["TARBP2_expression", "OS", "OS.time"]
)
df = df[df["OS"].isin([0, 1]) & (df["OS.time"] > 0)].copy()

median_expression = df["TARBP2_expression"].median()
df["High_Group"] = (
    df["TARBP2_expression"] > median_expression
).astype(int)

df["Expression_Group"] = np.where(
    df["High_Group"] == 1, "High TARBP2", "Low TARBP2"
)

# Age
age_col = "age_at_initial_pathologic_diagnosis"
if age_col in df.columns:
    df["Age"] = pd.to_numeric(df[age_col], errors="coerce")
else:
    df["Age"] = np.nan

# Gender: use Female as reference, Male = 1.
if "gender" in df.columns:
    gender = df["gender"].astype(str).str.strip().str.upper()
    df["Gender_Male"] = gender.map({
        "MALE": 1,
        "FEMALE": 0
    })
else:
    df["Gender_Male"] = np.nan

# Stage: prefer AJCC pathologic stage; use clinical stage if unavailable.
def simplify_stage(value):
    if pd.isna(value):
        return np.nan

    text = str(value).strip().upper()
    if text in {"", "NAN", "NA", "N/A", "NOT reported".upper(), "UNKNOWN"}:
        return np.nan

    # Detect stage I/II versus III/IV, including substages such as Stage IIA.
    match = re.search(r"STAGE\s*([IVX]+)", text)
    if not match:
        match = re.search(r"\b([IVX]+)[ABC]?\b", text)

    if not match:
        return np.nan

    roman = match.group(1)
    if roman in {"I", "II"}:
        return 0  # Early stage
    if roman in {"III", "IV"}:
        return 1  # Advanced stage
    return np.nan


# Prefer AJCC pathologic stage, but fill patient-level missing values
# with clinical stage when available.
df["Advanced_Stage"] = np.nan
stage_sources_used = []

if "ajcc_pathologic_tumor_stage" in df.columns:
    ajcc_stage = df["ajcc_pathologic_tumor_stage"].apply(simplify_stage)
    df["Advanced_Stage"] = ajcc_stage
    if ajcc_stage.notna().any():
        stage_sources_used.append("AJCC pathologic stage")

if "clinical_stage" in df.columns:
    clinical_stage = df["clinical_stage"].apply(simplify_stage)
    missing_stage = df["Advanced_Stage"].isna()
    df.loc[missing_stage, "Advanced_Stage"] = clinical_stage[missing_stage]
    if (missing_stage & clinical_stage.notna()).any():
        stage_sources_used.append("clinical stage（補足 AJCC 缺失）")

stage_source = " + ".join(stage_sources_used) if stage_sources_used else None


# =========================================================
# 3. Dataset overview
# =========================================================

st.header("1. KIRC Dataset Overview")

total_patients = len(df)
death_events = int(df["OS"].sum())
high_n = int((df["High_Group"] == 1).sum())
low_n = int((df["High_Group"] == 0).sum())

c1, c2, c3, c4 = st.columns(4)
c1.metric("Matched Patients", f"{total_patients:,}")
c2.metric("Death Events", f"{death_events:,}")
c3.metric("High TARBP2", f"{high_n:,}")
c4.metric("Low TARBP2", f"{low_n:,}")

st.write(f"TARBP2 median expression: **{median_expression:.3f}**")
st.caption(
    "高表現組定義為 TARBP2_expression 高於 KIRC 內部中位數；"
    "等於中位數者歸入低表現組。每位病患若有多筆 primary-tumor 樣本，"
    "先以其 TARBP2 表現量中位數代表該病患。"
)

st.divider()


# =========================================================
# 4. Expression distribution
# =========================================================

st.header("2. TARBP2 Expression Distribution")

fig, ax = plt.subplots(figsize=(7, 5))
groups = [
    df.loc[df["High_Group"] == 0, "TARBP2_expression"].dropna(),
    df.loc[df["High_Group"] == 1, "TARBP2_expression"].dropna()
]
ax.boxplot(
    groups,
    labels=["Low TARBP2", "High TARBP2"],
    showfliers=False
)
ax.set_ylabel("TARBP2 Expression")
ax.set_title("TARBP2 Expression in KIRC")
ax.grid(axis="y", alpha=0.25)
st.pyplot(fig)
plt.close(fig)

st.divider()


# =========================================================
# 5. Kaplan-Meier and Log-rank test
# =========================================================

st.header("3. Kaplan-Meier Survival Analysis")

high = df[df["High_Group"] == 1].copy()
low = df[df["High_Group"] == 0].copy()

kmf_high = KaplanMeierFitter()
kmf_low = KaplanMeierFitter()

fig, ax = plt.subplots(figsize=(9, 6))

kmf_high.fit(
    durations=high["OS.time"],
    event_observed=high["OS"],
    label="High TARBP2"
)
kmf_low.fit(
    durations=low["OS.time"],
    event_observed=low["OS"],
    label="Low TARBP2"
)

kmf_high.plot_survival_function(
    ax=ax, ci_show=True, show_censors=True
)
kmf_low.plot_survival_function(
    ax=ax, ci_show=True, show_censors=True
)

ax.set_title("Kaplan-Meier Overall Survival Curve - KIRC")
ax.set_xlabel("Time (Days)")
ax.set_ylabel("Overall Survival Probability")
ax.grid(alpha=0.25)
st.pyplot(fig)
plt.close(fig)

logrank_result = logrank_test(
    high["OS.time"],
    low["OS.time"],
    event_observed_A=high["OS"],
    event_observed_B=low["OS"]
)
logrank_p = float(logrank_result.p_value)

c1, c2, c3 = st.columns(3)
c1.metric("Log-rank P-value", f"{logrank_p:.4e}")
c2.metric("High-group deaths", f"{int(high['OS'].sum()):,}")
c3.metric("Low-group deaths", f"{int(low['OS'].sum()):,}")

if logrank_p < 0.05:
    st.success("Log-rank test：兩組 Overall Survival 曲線有統計顯著差異（P < 0.05）。")
else:
    st.info("Log-rank test：目前未觀察到兩組 Overall Survival 曲線有統計顯著差異（P ≥ 0.05）。")

st.divider()


# =========================================================
# 6. Univariate Cox regression
# =========================================================

st.header("4. Univariate Cox Regression")

univariate_data = df[["OS.time", "OS", "High_Group"]].dropna().copy()

try:
    cox_uni = CoxPHFitter()
    cox_uni.fit(
        univariate_data,
        duration_col="OS.time",
        event_col="OS"
    )

    uni_row = cox_uni.summary.loc["High_Group"]
    uni_hr = float(np.exp(uni_row["coef"]))
    uni_ci_low = float(np.exp(uni_row["coef lower 95%"]))
    uni_ci_high = float(np.exp(uni_row["coef upper 95%"]))
    uni_p = float(uni_row["p"])

    c1, c2, c3 = st.columns(3)
    c1.metric("Hazard Ratio", f"{uni_hr:.3f}")
    c2.metric("95% CI", f"{uni_ci_low:.3f} – {uni_ci_high:.3f}")
    c3.metric("Cox P-value", f"{uni_p:.4e}")

    st.dataframe(
        cox_uni.summary.reset_index(),
        use_container_width=True
    )
except Exception as e:
    st.warning(f"單變項 Cox regression 無法完成：{e}")
    uni_hr = uni_ci_low = uni_ci_high = uni_p = np.nan

st.divider()


# =========================================================
# 7. Multivariable Cox regression
# =========================================================

st.header("5. Multivariable Cox Regression")

st.markdown("""
此模型嘗試同時納入 TARBP2 高低表現組、診斷年齡、性別及癌症分期，
評估 TARBP2 與 OS 的關聯在控制這些臨床變項後是否仍存在。
""")

st.write(
    f"分期資料來源：**{stage_source if stage_source else '找不到可用分期欄位'}**"
)

covariates = [
    "OS.time", "OS", "High_Group",
    "Age", "Gender_Male", "Advanced_Stage"
]

multi_data = df[covariates].copy()
multi_data = multi_data.replace([np.inf, -np.inf], np.nan).dropna()

st.metric("Complete-case patients for multivariable model", f"{len(multi_data):,}")
st.caption(
    "多變項模型採 complete-case analysis：任一納入變項缺失的病患會被排除。"
)

multi_fit_ok = False
cox_multi = None

if (
    len(multi_data) >= 30
    and multi_data["OS"].sum() >= 10
    and multi_data["High_Group"].nunique() == 2
    and multi_data["Gender_Male"].nunique() == 2
    and multi_data["Advanced_Stage"].nunique() == 2
):
    try:
        cox_multi = CoxPHFitter()
        cox_multi.fit(
            multi_data,
            duration_col="OS.time",
            event_col="OS"
        )
        multi_fit_ok = True

        summary = cox_multi.summary.copy()
        summary["HR"] = np.exp(summary["coef"])
        summary["HR 95% CI Lower"] = np.exp(summary["coef lower 95%"])
        summary["HR 95% CI Upper"] = np.exp(summary["coef upper 95%"])

        display_cols = [
            "HR",
            "HR 95% CI Lower",
            "HR 95% CI Upper",
            "p"
        ]
        st.dataframe(
            summary[display_cols].rename(columns={"p": "P-value"}),
            use_container_width=True
        )

        if "High_Group" in summary.index:
            row = summary.loc["High_Group"]
            adj_hr = float(row["HR"])
            adj_low = float(row["HR 95% CI Lower"])
            adj_high = float(row["HR 95% CI Upper"])
            adj_p = float(row["p"])

            c1, c2, c3 = st.columns(3)
            c1.metric("Adjusted HR (High vs Low)", f"{adj_hr:.3f}")
            c2.metric("Adjusted 95% CI", f"{adj_low:.3f} – {adj_high:.3f}")
            c3.metric("Adjusted P-value", f"{adj_p:.4e}")

            if adj_p < 0.05 and (adj_low > 1 or adj_high < 1):
                st.success(
                    "控制年齡、性別及分期後，TARBP2 高低表現與 OS 的關聯仍達統計顯著。"
                )
            else:
                st.info(
                    "控制年齡、性別及分期後，TARBP2 高低表現與 OS 的關聯未達統計顯著。"
                )

    except Exception as e:
        st.warning(f"多變項 Cox regression 無法完成：{e}")
else:
    st.warning(
        "目前完整臨床資料不足，或某個臨床變項沒有足夠變異，"
        "因此暫不執行多變項 Cox。請先檢查下方缺失值與變項分布。"
    )

with st.expander("查看臨床變項可用資料數"):
    availability = pd.DataFrame({
        "Variable": ["Age", "Gender", "Stage"],
        "Available N": [
            int(df["Age"].notna().sum()),
            int(df["Gender_Male"].notna().sum()),
            int(df["Advanced_Stage"].notna().sum())
        ],
        "Missing N": [
            int(df["Age"].isna().sum()),
            int(df["Gender_Male"].isna().sum()),
            int(df["Advanced_Stage"].isna().sum())
        ]
    })
    st.dataframe(availability, use_container_width=True)

st.divider()


# =========================================================
# 8. Proportional hazards assumption
# =========================================================

st.header("6. Proportional Hazards Assumption")

if multi_fit_ok and cox_multi is not None:
    try:
        ph_test = proportional_hazard_test(
            cox_multi,
            multi_data,
            time_transform="rank"
        )
        ph_result = ph_test.summary.reset_index()
        ph_result = ph_result.rename(
            columns={
                "index": "Variable",
                "test_statistic": "Test Statistic",
                "p": "P-value"
            }
        )
        st.dataframe(ph_result, use_container_width=True)

        if (ph_test.summary["p"] < 0.05).any():
            st.warning(
                "至少一個變項的比例風險假設檢定 P < 0.05，"
                "表示可能違反 proportional hazards assumption；"
                "需進一步檢查模型或考慮分層/時間變動效果。"
            )
        else:
            st.success(
                "目前檢定未發現明確違反比例風險假設的證據（所有 P ≥ 0.05）。"
            )
    except Exception as e:
        st.warning(f"比例風險假設檢查失敗：{e}")
else:
    st.info("多變項 Cox 模型未成功建立，因此暫時無法檢查比例風險假設。")


# =========================================================
# 9. Stratified Cox regression (stratified by cancer stage)
# =========================================================

st.header("7. Stratified Cox Regression")

st.markdown("""
前面的多變項 Cox 模型顯示 Advanced_Stage 可能不符合比例風險假設。
因此，這裡將 **Advanced_Stage 作為分層變項（strata）**，
允許早期與晚期病患具有不同的 baseline hazard function，
並重新估計 TARBP2、年齡及性別與 OS 的關聯。

注意：分層 Cox 不會估計 Advanced_Stage 本身的 HR；
此模型也不代表已經解決所有可能的偏差或混雜。
""")

stratified_fit_ok = False
cox_stratified = None

if (
    len(multi_data) >= 30
    and multi_data["OS"].sum() >= 10
    and multi_data["High_Group"].nunique() == 2
    and multi_data["Gender_Male"].nunique() == 2
    and multi_data["Advanced_Stage"].nunique() == 2
):
    try:
        stratified_data = multi_data[
            ["OS.time", "OS", "High_Group", "Age",
             "Gender_Male", "Advanced_Stage"]
        ].copy()

        cox_stratified = CoxPHFitter()
        cox_stratified.fit(
            stratified_data,
            duration_col="OS.time",
            event_col="OS",
            strata=["Advanced_Stage"]
        )
        stratified_fit_ok = True

        strat_summary = cox_stratified.summary.copy()
        strat_summary["HR"] = np.exp(strat_summary["coef"])
        strat_summary["HR 95% CI Lower"] = np.exp(
            strat_summary["coef lower 95%"]
        )
        strat_summary["HR 95% CI Upper"] = np.exp(
            strat_summary["coef upper 95%"]
        )

        st.write("分層變項：**Advanced_Stage（早期／晚期）**")
        st.write(
            f"分層 Cox 使用完整案例數：**{len(stratified_data):,}**；"
            f"死亡事件數：**{int(stratified_data['OS'].sum()):,}**"
        )

        strat_display_cols = [
            "HR", "HR 95% CI Lower", "HR 95% CI Upper", "p"
        ]
        st.dataframe(
            strat_summary[strat_display_cols].rename(
                columns={"p": "P-value"}
            ),
            use_container_width=True
        )

        if "High_Group" in strat_summary.index:
            strat_row = strat_summary.loc["High_Group"]
            strat_hr = float(strat_row["HR"])
            strat_low = float(strat_row["HR 95% CI Lower"])
            strat_high = float(strat_row["HR 95% CI Upper"])
            strat_p = float(strat_row["p"])

            c1, c2, c3 = st.columns(3)
            c1.metric("Stratified HR (High vs Low)", f"{strat_hr:.3f}")
            c2.metric(
                "Stratified 95% CI",
                f"{strat_low:.3f} – {strat_high:.3f}"
            )
            c3.metric("Stratified P-value", f"{strat_p:.4e}")

            if strat_p < 0.05 and (strat_low > 1 or strat_high < 1):
                st.success(
                    "以癌症分期分層後，TARBP2 高低表現與 OS 的關聯仍達統計顯著。"
                )
            else:
                st.info(
                    "以癌症分期分層後，TARBP2 高低表現與 OS 的關聯未達統計顯著。"
                )

        st.subheader("Stratified Cox：比例風險假設檢查")
        strat_ph_test = proportional_hazard_test(
            cox_stratified,
            stratified_data,
            time_transform="rank"
        )
        strat_ph_result = strat_ph_test.summary.reset_index()
        strat_ph_result = strat_ph_result.rename(
            columns={
                "index": "Variable",
                "test_statistic": "Test Statistic",
                "p": "P-value"
            }
        )
        st.dataframe(strat_ph_result, use_container_width=True)

        if (strat_ph_test.summary["p"] < 0.05).any():
            violating = strat_ph_test.summary.index[
                strat_ph_test.summary["p"] < 0.05
            ].tolist()
            st.warning(
                "分層 Cox 中仍有變項可能違反比例風險假設："
                + "、".join(map(str, violating))
                + "。需審慎解讀，必要時考慮時間變動效果或其他模型。"
            )
        else:
            st.success(
                "分層 Cox 的 PH 檢定未發現明確違反比例風險假設的證據"
                "（所有 P ≥ 0.05）。"
            )

        st.caption(
            "分層 Cox 的 HR 是在不同 Advanced_Stage 分層內估計的共同協變項效果；"
            "Advanced_Stage 本身不會出現在 HR 表格中。"
        )

    except Exception as e:
        st.warning(f"分層 Cox regression 或其 PH 檢查無法完成：{e}")
else:
    st.info(
        "目前完整臨床資料不足，或分層／協變項沒有足夠變異，"
        "因此無法執行分層 Cox regression。"
    )

st.divider()


# =========================================================
# 10. Interpretation notes
# =========================================================


st.header("8. Interpretation Notes")

st.markdown("""
- **HR > 1**：高 TARBP2 表現組的估計事件風險高於低表現組。
- **HR < 1**：高 TARBP2 表現組的估計事件風險低於低表現組。
- **95% CI 包含 1**：通常表示 HR 與 1 無差異的情況仍與資料相容。
- **P-value < 0.05**：常用統計判讀門檻，不等於生物學重要性。
- 多變項 Cox 的結果取決於納入變項、缺失資料處理及模型假設。
- 本分析為觀察性資料分析，不能據此推論 TARBP2 對存活具有因果作用。
- 本模組僅供研究與教學用途，不可用於個人臨床預後判斷。
""")

st.caption(
    "資料來源：TCGA Pan-Cancer TARBP2 primary-tumor expression "
    "及 Survival_SupplementalTable_S1_20171025_xena_sp。"
)
