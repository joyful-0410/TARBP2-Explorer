import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test


# =========================================================
# TARBP2 Pan-Cancer Survival Analysis
# =========================================================

st.set_page_config(
    page_title="TARBP2 Pan-Cancer Survival Analysis",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 TARBP2 Pan-Cancer Survival Analysis")

st.markdown("""
本模組使用 TCGA Pan-Cancer TARBP2 RNA-seq 表現量與臨床存活資料，
逐一分析不同癌症中 TARBP2 高、低表現組的 Overall Survival（OS）。
每一種癌症分別進行 Log-rank test 與 Cox proportional hazards regression，
並使用 Benjamini-Hochberg 方法進行多重比較校正。
""")

st.divider()


# =========================================================
# 1. Load and merge data
# =========================================================

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
        subset=["sample", "TARBP2_expression", "_primary_disease"]
    )
    clinical = clinical.dropna(
        subset=[
            "sample",
            "_PATIENT",
            "cancer type abbreviation",
            "OS",
            "OS.time"
        ]
    )

    clinical = clinical[
        clinical["OS"].isin([0, 1]) &
        (clinical["OS.time"] > 0)
    ].copy()

    expression = expression.drop_duplicates(
        subset=["sample"], keep="first"
    )
    clinical = clinical.drop_duplicates(
        subset=["sample"], keep="first"
    )

    merged = expression.merge(
        clinical,
        on="sample",
        how="inner",
        validate="one_to_one"
    )

    # 本分析以病患為單位；若同一病患有多筆樣本，只保留一筆。
    merged = merged.sort_values("sample")
    merged = merged.drop_duplicates(
        subset=["_PATIENT"], keep="first"
    )

    if merged.empty:
        raise ValueError("Expression 與 Clinical 資料沒有成功配對。")

    return merged


try:
    df = load_survival_data()
except Exception as e:
    st.error(f"資料讀取或合併失敗：{e}")
    st.stop()


# =========================================================
# 2. Benjamini-Hochberg FDR correction
# =========================================================

def benjamini_hochberg(p_values):
    """回傳與輸入順序相同的 Benjamini-Hochberg adjusted p-values。"""
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full(p_values.shape, np.nan, dtype=float)

    valid_mask = np.isfinite(p_values)
    valid_p = p_values[valid_mask]

    if len(valid_p) == 0:
        return adjusted

    n = len(valid_p)
    order = np.argsort(valid_p)
    ranked_p = valid_p[order]

    ranked_adjusted = ranked_p * n / np.arange(1, n + 1)
    ranked_adjusted = np.minimum.accumulate(ranked_adjusted[::-1])[::-1]
    ranked_adjusted = np.clip(ranked_adjusted, 0, 1)

    restored = np.empty(n, dtype=float)
    restored[order] = ranked_adjusted
    adjusted[valid_mask] = restored

    return adjusted


# =========================================================
# 3. Analyze each cancer type
# =========================================================

def analyze_one_cancer(cancer_df, cancer_code, min_patients):
    cancer_df = cancer_df.copy()

    cancer_df["TARBP2_expression"] = pd.to_numeric(
        cancer_df["TARBP2_expression"], errors="coerce"
    )
    cancer_df["OS"] = pd.to_numeric(cancer_df["OS"], errors="coerce")
    cancer_df["OS.time"] = pd.to_numeric(
        cancer_df["OS.time"], errors="coerce"
    )

    cancer_df = cancer_df.dropna(
        subset=["TARBP2_expression", "OS", "OS.time"]
    )
    cancer_df = cancer_df[
        cancer_df["OS"].isin([0, 1]) & (cancer_df["OS.time"] > 0)
    ].copy()

    n_total = len(cancer_df)
    if n_total < min_patients:
        return None

    median_expression = cancer_df["TARBP2_expression"].median()

    # 高表現組定義為高於該癌症內部的中位數；其餘歸入低表現組。
    cancer_df["High_Group"] = (
        cancer_df["TARBP2_expression"] > median_expression
    ).astype(int)

    high = cancer_df[cancer_df["High_Group"] == 1].copy()
    low = cancer_df[cancer_df["High_Group"] == 0].copy()

    if len(high) < 2 or len(low) < 2:
        return None

    high_events = int(high["OS"].sum())
    low_events = int(low["OS"].sum())

    # Log-rank test
    try:
        logrank_result = logrank_test(
            high["OS.time"],
            low["OS.time"],
            event_observed_A=high["OS"],
            event_observed_B=low["OS"]
        )
        logrank_p = float(logrank_result.p_value)
    except Exception:
        logrank_p = np.nan

    # Cox model: High_Group = 1 means high TARBP2 expression.
    hr = np.nan
    ci_lower = np.nan
    ci_upper = np.nan
    cox_p = np.nan

    cox_data = cancer_df[["OS.time", "OS", "High_Group"]].copy()

    try:
        # Cox regression requires sufficient event information in the dataset.
        if int(cox_data["OS"].sum()) >= 2:
            cph = CoxPHFitter()
            cph.fit(
                cox_data,
                duration_col="OS.time",
                event_col="OS"
            )

            hr = float(cph.hazard_ratios_["High_Group"])
            ci_lower = float(
                cph.confidence_intervals_.loc["High_Group"].iloc[0]
            )
            ci_upper = float(
                cph.confidence_intervals_.loc["High_Group"].iloc[1]
            )
            cox_p = float(cph.summary.loc["High_Group", "p"])

            # lifelines confidence_intervals_ is on the log(HR) scale.
            ci_lower = float(np.exp(ci_lower))
            ci_upper = float(np.exp(ci_upper))

    except Exception:
        pass

    return {
        "Cancer": cancer_code,
        "Patients": n_total,
        "High N": len(high),
        "Low N": len(low),
        "High Deaths": high_events,
        "Low Deaths": low_events,
        "Median Expression": median_expression,
        "Log-rank P": logrank_p,
        "HR": hr,
        "HR 95% CI Lower": ci_lower,
        "HR 95% CI Upper": ci_upper,
        "Cox P": cox_p
    }


st.subheader("1. Dataset Overview")

col1, col2, col3 = st.columns(3)
col1.metric("Matched Patients", f"{len(df):,}")
col2.metric("Death Events", f"{int(df['OS'].sum()):,}")
col3.metric(
    "Cancer Types",
    int(df["cancer type abbreviation"].nunique())
)

st.caption(
    "分析單位為病患；同一病患若有多筆樣本，只保留一筆。"
    "OS：1 = death event，0 = censored/alive。"
)

st.divider()


# =========================================================
# 4. Analysis settings and run analysis
# =========================================================

st.subheader("2. Pan-Cancer Analysis")

min_patients = st.number_input(
    "每種癌症最低病患數",
    min_value=10,
    max_value=200,
    value=20,
    step=5,
    help="低於此病患數的癌症不納入分析。"
)

run_analysis = st.button(
    "開始分析全部癌症",
    type="primary",
    use_container_width=True
)

if "pan_cancer_results" not in st.session_state:
    st.session_state["pan_cancer_results"] = None

if run_analysis:
    result_rows = []
    cancer_codes = sorted(
        df["cancer type abbreviation"].dropna().unique().tolist()
    )

    progress = st.progress(0)
    status = st.empty()

    for i, cancer_code in enumerate(cancer_codes):
        status.write(
            f"正在分析 {cancer_code}（{i + 1}/{len(cancer_codes)}）..."
        )

        subset = df[
            df["cancer type abbreviation"] == cancer_code
        ].copy()

        result = analyze_one_cancer(
            subset,
            cancer_code,
            int(min_patients)
        )

        if result is not None:
            result_rows.append(result)

        progress.progress((i + 1) / len(cancer_codes))

    progress.empty()
    status.empty()

    if not result_rows:
        st.session_state["pan_cancer_results"] = None
        st.error(
            "沒有癌症符合分析條件。請降低最低病患數，"
            "或檢查臨床資料與表現量資料。"
        )
    else:
        results = pd.DataFrame(result_rows)

        results["Log-rank FDR"] = benjamini_hochberg(
            results["Log-rank P"].values
        )
        results["Cox FDR"] = benjamini_hochberg(
            results["Cox P"].values
        )

        results = results.sort_values(
            "Log-rank P",
            na_position="last"
        ).reset_index(drop=True)

        st.session_state["pan_cancer_results"] = results


# =========================================================
# 5. Display results
# =========================================================

results = st.session_state["pan_cancer_results"]

if results is not None:
    st.divider()
    st.subheader("3. Analysis Results")

    analyzed_count = len(results)
    logrank_significant = int(
        (results["Log-rank FDR"] < 0.05).sum()
    )
    cox_significant = int(
        (results["Cox FDR"] < 0.05).sum()
    )

    col1, col2, col3 = st.columns(3)
    col1.metric("Analyzed Cancer Types", analyzed_count)
    col2.metric("Log-rank FDR < 0.05", logrank_significant)
    col3.metric("Cox FDR < 0.05", cox_significant)

    st.caption(
        "FDR 分別針對所有可計算的 Log-rank P-values 與 Cox P-values "
        "進行 Benjamini-Hochberg 校正。"
    )

    display_columns = [
        "Cancer",
        "Patients",
        "High N",
        "Low N",
        "High Deaths",
        "Low Deaths",
        "HR",
        "HR 95% CI Lower",
        "HR 95% CI Upper",
        "Log-rank P",
        "Log-rank FDR",
        "Cox P",
        "Cox FDR"
    ]

    st.dataframe(
        results[display_columns].style.format({
            "HR": "{:.3f}",
            "HR 95% CI Lower": "{:.3f}",
            "HR 95% CI Upper": "{:.3f}",
            "Log-rank P": "{:.4g}",
            "Log-rank FDR": "{:.4g}",
            "Cox P": "{:.4g}",
            "Cox FDR": "{:.4g}"
        }, na_rep="NA"),
        use_container_width=True,
        hide_index=True
    )

    csv_data = results.to_csv(index=False).encode("utf-8-sig")
    st.download_button(
        "下載 Pan-Cancer Survival 分析結果 CSV",
        data=csv_data,
        file_name="TARBP2_PanCancer_Survival_Results.csv",
        mime="text/csv"
    )

    st.divider()


    # =====================================================
    # 6. Forest plot
    # =====================================================

    st.subheader("4. Cox Regression Forest Plot")

    forest_df = results.dropna(
        subset=["HR", "HR 95% CI Lower", "HR 95% CI Upper"]
    ).copy()

    forest_df = forest_df[
        (forest_df["HR"] > 0) &
        (forest_df["HR 95% CI Lower"] > 0) &
        (forest_df["HR 95% CI Upper"] > 0)
    ].sort_values("HR")

    if forest_df.empty:
        st.warning(
            "目前沒有可繪製的 Cox regression 結果。"
            "可能是部分癌症的 Cox model 無法穩定估計。"
        )
    else:
        fig_height = max(5, 0.42 * len(forest_df) + 1.5)
        fig, ax = plt.subplots(figsize=(10, fig_height))

        y_positions = np.arange(len(forest_df))
        hr_values = forest_df["HR"].to_numpy()
        lower_errors = (
            hr_values - forest_df["HR 95% CI Lower"].to_numpy()
        )
        upper_errors = (
            forest_df["HR 95% CI Upper"].to_numpy() - hr_values
        )

        ax.errorbar(
            hr_values,
            y_positions,
            xerr=np.vstack([lower_errors, upper_errors]),
            fmt="o",
            capsize=3
        )

        ax.axvline(
            x=1,
            linestyle="--",
            linewidth=1
        )

        ax.set_yticks(y_positions)
        ax.set_yticklabels(forest_df["Cancer"])
        ax.set_xlabel("Hazard Ratio (High vs Low TARBP2 Expression)")
        ax.set_ylabel("Cancer Type")
        ax.set_title("TARBP2 Pan-Cancer Cox Regression")
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()

        st.pyplot(fig)
        plt.close(fig)

        st.caption(
            "HR > 1：高 TARBP2 表現組的估計事件風險較高；"
            "HR < 1：估計事件風險較低。"
            "若 95% CI 跨過 1，通常表示未達統計顯著。"
        )

    st.divider()


    # =====================================================
    # 7. Significant findings
    # =====================================================

    st.subheader("5. Significant Findings")

    significant = results[
        results["Log-rank FDR"] < 0.05
    ].copy()

    if significant.empty:
        st.info(
            "目前沒有癌症在 Log-rank test 經 FDR 多重比較校正後 "
            "達到 0.05 的顯著門檻。"
        )
    else:
        st.success(
            f"共有 {len(significant)} 種癌症在 Log-rank test "
            "經 FDR 校正後達到顯著。"
        )
        st.dataframe(
            significant[
                [
                    "Cancer",
                    "Patients",
                    "Log-rank P",
                    "Log-rank FDR",
                    "HR",
                    "HR 95% CI Lower",
                    "HR 95% CI Upper",
                    "Cox P",
                    "Cox FDR"
                ]
            ],
            use_container_width=True,
            hide_index=True
        )

    st.divider()

    st.subheader("6. Interpretation Notes")
    st.markdown("""
    - HR（Hazard Ratio）比較高 TARBP2 表現組相對於低表現組的估計事件風險。
    - HR > 1 表示高表現組估計風險較高；HR < 1 表示估計風險較低。
    - 95% 信賴區間若包含 1，代表該 HR 估計的不確定性仍涵蓋無差異的情況。
    - P-value 與 FDR < 0.05 是常用統計判讀門檻，但不能單獨代表生物學重要性。
    - 本分析為觀察性資料分析，不能據此推論 TARBP2 對存活結果具有因果作用。
    - 本模組僅供研究與教學用途，不可用於個人臨床預後判斷。
    """)

else:
    st.info(
        "按下「開始分析全部癌症」後，系統會逐一執行各癌症的 "
        "Log-rank test 與 Cox regression，並產生結果表格及 Forest Plot。"
    )
