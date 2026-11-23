import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

try:
    from scipy.stats import mannwhitneyu, wilcoxon
except ImportError:
    mannwhitneyu = None
    wilcoxon = None

try:
    from lifelines import CoxPHFitter
    from lifelines.statistics import logrank_test
except ImportError:
    CoxPHFitter = None
    logrank_test = None


# =========================================================
# TARBP2 Role Assessment
# Evidence integration: expression, tumor-normal, survival, STRING
# =========================================================

st.set_page_config(
    page_title="TARBP2 Role Assessment",
    page_icon="🧭",
    layout="wide"
)

st.title("🧭 TARBP2 Role Assessment")
st.markdown("""
本頁整合 TCGA 表現量、腫瘤與正常組織比較、生存分析及 STRING 網絡資料，
整理 TARBP2 在癌症中的候選角色證據。本頁屬於**證據整合與假說產生**，
不代表已證實 TARBP2 具有促癌或抑癌的因果作用。
""")

st.info(
    "判讀原則：表現量差異不等於功能改變；生存關聯不等於因果；"
    "STRING functional association 也不等於直接物理結合。"
)

BASE = Path(".")
CLINICAL_PATH = BASE / "data" / "clinical" / "Survival_SupplementalTable_S1_20171025_xena_sp"


def find_existing(candidates):
    for item in candidates:
        path = Path(item)
        if path.exists() and path.is_file():
            return path
    return None


def read_csv_safe(path):
    if path is None:
        return None
    try:
        return pd.read_csv(path)
    except Exception as exc:
        st.warning(f"無法讀取 {path}: {exc}")
        return None


def patient_id_from_barcode(series):
    """TCGA barcode 前 12 個字元作為病患層級識別碼。"""
    ids = series.astype(str).str.strip().str[:12]
    return ids.where(~ids.isin(["", "nan", "None"]), np.nan)


def benjamini_hochberg(p_values):
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    p = values[valid]
    if len(p) == 0:
        return adjusted
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    restored = np.empty(len(p), dtype=float)
    restored[order] = np.clip(q, 0, 1)
    adjusted[valid] = restored
    return adjusted


def load_survival():
    if CoxPHFitter is None or logrank_test is None:
        raise ImportError("缺少 lifelines，請先安裝：pip install lifelines")

    expression_path = find_existing([
        "data/TCGA/TARBP2_primary_tumor.csv",
        "data/expression/TARBP2_primary_tumor.csv"
    ])
    if expression_path is None or not CLINICAL_PATH.exists():
        raise FileNotFoundError("找不到 TARBP2_primary_tumor.csv 或臨床生存資料。")

    exp = pd.read_csv(expression_path)
    cli = pd.read_csv(CLINICAL_PATH, sep="\t")

    exp_need = ["sample", "TARBP2_expression", "_primary_disease"]
    cli_need = ["sample", "_PATIENT", "cancer type abbreviation", "OS", "OS.time"]

    for col in exp_need:
        if col not in exp.columns:
            raise ValueError(f"表現量資料缺少欄位：{col}")
    for col in cli_need:
        if col not in cli.columns:
            raise ValueError(f"臨床資料缺少欄位：{col}")

    exp = exp.copy()
    cli = cli.copy()
    exp["sample"] = exp["sample"].astype(str).str.strip()
    cli["sample"] = cli["sample"].astype(str).str.strip()
    cli["_PATIENT"] = cli["_PATIENT"].astype(str).str.strip()
    exp["TARBP2_expression"] = pd.to_numeric(
        exp["TARBP2_expression"], errors="coerce"
    )
    cli["OS"] = pd.to_numeric(cli["OS"], errors="coerce")
    cli["OS.time"] = pd.to_numeric(cli["OS.time"], errors="coerce")

    exp = exp.dropna(subset=["sample", "TARBP2_expression", "_primary_disease"])
    cli = cli.dropna(subset=cli_need)
    cli = cli[cli["OS"].isin([0, 1]) & (cli["OS.time"] > 0)].copy()

    # Avoid accidental many-to-many merge at sample level.
    exp = exp.drop_duplicates("sample", keep="first")
    cli = cli.drop_duplicates("sample", keep="first")

    merged = exp.merge(
        cli, on="sample", how="inner", validate="one_to_one"
    )
    merged["_PATIENT"] = merged["_PATIENT"].astype(str).str.strip()
    merged = merged[~merged["_PATIENT"].isin(["", "nan", "None"])].copy()

    if merged.empty:
        raise ValueError("表現量與臨床資料合併後沒有可用樣本。")

    # One record per patient: median expression; retain one clinical record.
    # Clinical OS fields should be patient-level and are checked for consistency.
    clinical_cols = [
        c for c in [
            "_PATIENT", "cancer type abbreviation", "OS", "OS.time",
            "age_at_initial_pathologic_diagnosis", "gender",
            "ajcc_pathologic_tumor_stage", "clinical_stage"
        ] if c in merged.columns
    ]
    clinical = merged[clinical_cols].drop_duplicates("_PATIENT", keep="first")
    expression_patient = (
        merged.groupby("_PATIENT", as_index=False)
        .agg(TARBP2_expression=("TARBP2_expression", "median"))
    )
    patient_df = expression_patient.merge(
        clinical, on="_PATIENT", how="inner", validate="one_to_one"
    )

    if patient_df["_PATIENT"].duplicated().any():
        raise ValueError("病患層級整理後仍有重複 _PATIENT。")

    return patient_df


def tumor_normal_by_cancer(tumor_df, normal_df):
    """Cancer-wise comparison; paired Wilcoxon when >=5 patient pairs."""
    required = {"TARBP2_expression", "_primary_disease", "sample_type"}
    if not required.issubset(tumor_df.columns):
        raise ValueError(
            "腫瘤資料缺少欄位：" + ", ".join(sorted(required - set(tumor_df.columns)))
        )
    if not required.issubset(normal_df.columns):
        raise ValueError(
            "正常資料缺少欄位：" + ", ".join(sorted(required - set(normal_df.columns)))
        )

    tumor = tumor_df.copy()
    normal = normal_df.copy()
    for frame in [tumor, normal]:
        frame["TARBP2_expression"] = pd.to_numeric(
            frame["TARBP2_expression"], errors="coerce"
        )
        frame["_primary_disease"] = frame["_primary_disease"].astype(str)

    tumor = tumor[tumor["sample_type"].astype(str).str.contains(
        "Primary Tumor", case=False, na=False
    )]
    normal = normal[normal["sample_type"].astype(str).str.contains(
        "Solid Tissue Normal", case=False, na=False
    )]

    has_barcode = "sample" in tumor.columns and "sample" in normal.columns
    if has_barcode:
        tumor["Patient_ID"] = patient_id_from_barcode(tumor["sample"])
        normal["Patient_ID"] = patient_id_from_barcode(normal["sample"])

    rows = []
    diseases = sorted(
        set(tumor["_primary_disease"].dropna())
        .intersection(set(normal["_primary_disease"].dropna()))
    )

    for disease in diseases:
        t = tumor[tumor["_primary_disease"] == disease].dropna(
            subset=["TARBP2_expression"]
        ).copy()
        n = normal[normal["_primary_disease"] == disease].dropna(
            subset=["TARBP2_expression"]
        ).copy()

        if len(t) < 5 or len(n) < 5:
            continue

        pair_count = 0
        method = "Mann–Whitney U（獨立樣本）"
        note = "配對數不足或無法辨識病患 ID，採獨立樣本比較。"

        if has_barcode:
            t_patient = (
                t.dropna(subset=["Patient_ID"])
                .groupby("Patient_ID")["TARBP2_expression"].median()
            )
            n_patient = (
                n.dropna(subset=["Patient_ID"])
                .groupby("Patient_ID")["TARBP2_expression"].median()
            )
            common = t_patient.index.intersection(n_patient.index)
            pair_count = len(common)

        if has_barcode and pair_count >= 5:
            tv = t_patient.loc[common]
            nv = n_patient.loc[common]
            differences = tv - nv
            if wilcoxon is None:
                continue
            if np.allclose(differences.to_numpy(), 0):
                statistic, p_value = 0.0, 1.0
            else:
                statistic, p_value = wilcoxon(tv, nv, alternative="two-sided")
            tumor_median = tv.median()
            normal_median = nv.median()
            tumor_used = normal_used = pair_count
            method = "Wilcoxon signed-rank（配對樣本）"
            note = "同一病患腫瘤與正常樣本配對；同類多筆樣本以中位數彙整。"
        else:
            if mannwhitneyu is None:
                continue
            tv = t["TARBP2_expression"]
            nv = n["TARBP2_expression"]
            statistic, p_value = mannwhitneyu(
                tv, nv, alternative="two-sided"
            )
            tumor_median = tv.median()
            normal_median = nv.median()
            tumor_used, normal_used = len(tv), len(nv)

        rows.append({
            "Cancer Type": disease,
            "Tumor N": len(t),
            "Normal N": len(n),
            "Matched Pairs": pair_count,
            "Tumor N Used": tumor_used,
            "Normal N Used": normal_used,
            "Tumor Median": tumor_median,
            "Normal Median": normal_median,
            "Median Difference": tumor_median - normal_median,
            "Test Method": method,
            "P-value": float(p_value),
            "Analysis Note": note
        })

    result = pd.DataFrame(rows)
    if not result.empty:
        result["FDR"] = benjamini_hochberg(result["P-value"])
        result = result.sort_values("FDR")
    return result


# ---------------------------------------------------------
# 1. Expression evidence
# ---------------------------------------------------------
st.header("1. Pan-Cancer Expression Evidence")

expression_path = find_existing([
    "data/TCGA/TARBP2_primary_tumor.csv",
    "data/expression/TARBP2_primary_tumor.csv"
])
summary_path = find_existing([
    "data/TCGA/TARBP2_cancer_summary.csv",
    "data/expression/TARBP2_cancer_summary.csv"
])
normal_path = find_existing([
    "data/TCGA/TARBP2_normal_tissue.csv",
    "data/expression/TARBP2_normal_tissue.csv"
])

expr = read_csv_safe(expression_path)
summary = read_csv_safe(summary_path)
normal = read_csv_safe(normal_path)

if expr is None:
    st.warning("找不到 TCGA TARBP2 primary tumor 表現量資料。請確認資料路徑。")
else:
    required = {"TARBP2_expression", "_primary_disease"}
    if not required.issubset(expr.columns):
        st.error(f"表現量資料缺少欄位：{sorted(required - set(expr.columns))}")
        expr = None
    else:
        expr["TARBP2_expression"] = pd.to_numeric(
            expr["TARBP2_expression"], errors="coerce"
        )
        expr = expr.dropna(subset=["TARBP2_expression", "_primary_disease"])
        c1, c2, c3 = st.columns(3)
        c1.metric("Primary tumor samples", f"{len(expr):,}")
        c2.metric("Cancer types", int(expr["_primary_disease"].nunique()))
        c3.metric(
            "Overall median expression",
            f"{expr['TARBP2_expression'].median():.3f}"
        )

        if summary is not None and {
            "_primary_disease", "median"
        }.issubset(summary.columns):
            show_cols = [
                c for c in [
                    "_primary_disease", "sample_count", "mean", "median", "std"
                ] if c in summary.columns
            ]
            display_summary = summary[show_cols].sort_values(
                "median", ascending=False
            )
        else:
            display_summary = (
                expr.groupby("_primary_disease")["TARBP2_expression"]
                .agg(
                    sample_count="count", mean="mean",
                    median="median", std="std"
                )
                .reset_index()
                .sort_values("median", ascending=False)
            )

        st.subheader("癌症別 TARBP2 表現量摘要")
        st.dataframe(
            display_summary, use_container_width=True, hide_index=True
        )
        st.caption(
            "此處為各癌症 primary tumor 樣本的描述性統計。"
            "不同癌症的樣本組成與資料分布可能不同，不能僅以跨癌症中位數高低推論癌症風險。"
        )


# ---------------------------------------------------------
# 2. Tumor-normal evidence
# ---------------------------------------------------------
st.divider()
st.header("2. Tumor vs Normal Evidence")

if expr is None:
    st.info("需先成功載入 primary tumor 表現量資料，才能進行腫瘤／正常比較。")
elif normal is None:
    st.warning(
        "找不到 TARBP2_normal_tissue.csv，因此本區無法比較。"
        "缺少正常組資料不代表沒有表現差異。"
    )
else:
    try:
        tn_results = tumor_normal_by_cancer(expr, normal)
        if tn_results.empty:
            st.warning(
                "目前沒有癌症種類同時具備足夠腫瘤與正常樣本，"
                "或缺少必要欄位／統計套件。"
            )
        else:
            paired_count = int(
                tn_results["Test Method"].str.contains("配對樣本").sum()
            )
            sig_count = int((tn_results["FDR"] < 0.05).sum())
            c1, c2, c3 = st.columns(3)
            c1.metric("Cancer types analyzed", len(tn_results))
            c2.metric("Cancer types with paired test", paired_count)
            c3.metric("FDR < 0.05", sig_count)

            st.dataframe(
                tn_results.style.format({
                    "Tumor Median": "{:.3f}",
                    "Normal Median": "{:.3f}",
                    "Median Difference": "{:.3f}",
                    "P-value": "{:.4g}",
                    "FDR": "{:.4g}"
                }),
                use_container_width=True,
                hide_index=True
            )
            st.download_button(
                "下載癌症別 Tumor vs Normal 證據 CSV",
                tn_results.to_csv(index=False).encode("utf-8-sig"),
                file_name="TARBP2_Role_Assessment_Tumor_Normal.csv",
                mime="text/csv"
            )
            st.caption(
                "若同一病患同時有腫瘤與正常樣本且配對數至少 5，使用 Wilcoxon signed-rank；"
                "否則使用 Mann–Whitney U。各癌症 P-value 以 BH 方法進行 FDR 校正。"
                "此分析未校正臨床共變項。"
            )
    except Exception as exc:
        st.error(f"Tumor vs Normal 分析失敗：{exc}")


# ---------------------------------------------------------
# 3. Survival evidence
# ---------------------------------------------------------
st.divider()
st.header("3. Prognostic Evidence")

st.markdown(
    "以每種癌症內 TARBP2 表現量中位數分組，進行 Log-rank 與單變項 Cox 分析。"
    "這些結果是未校正臨床共變項的探索性關聯，不應直接稱為獨立預後因子。"
)

try:
    survival_df = load_survival()
    st.caption(
        f"成功配對 {len(survival_df):,} 位病患；"
        f"死亡事件 {int(survival_df['OS'].sum()):,}；"
        f"癌症種類 {survival_df['cancer type abbreviation'].nunique()}。"
    )

    min_n = st.number_input(
        "每種癌症最低病患數",
        min_value=10, max_value=200, value=20, step=5,
        key="role_min_n"
    )

    if st.button(
        "重新計算 Pan-Cancer Survival Evidence",
        type="primary"
    ):
        rows = []
        errors = []

        for cancer in sorted(
            survival_df["cancer type abbreviation"].dropna().unique()
        ):
            sub = survival_df[
                survival_df["cancer type abbreviation"] == cancer
            ].copy()
            sub = sub.dropna(
                subset=["TARBP2_expression", "OS", "OS.time"]
            )

            if len(sub) < int(min_n):
                errors.append({
                    "Cancer": cancer,
                    "Reason": f"病患數低於門檻（n={len(sub)}）"
                })
                continue

            median_value = sub["TARBP2_expression"].median()
            sub["High_Group"] = (
                sub["TARBP2_expression"] > median_value
            ).astype(int)
            high = sub[sub["High_Group"] == 1]
            low = sub[sub["High_Group"] == 0]

            if (
                len(high) < 2 or len(low) < 2
                or int(sub["OS"].sum()) < 2
            ):
                errors.append({
                    "Cancer": cancer,
                    "Reason": "高低組樣本數或死亡事件不足"
                })
                continue

            try:
                lr = logrank_test(
                    high["OS.time"], low["OS.time"],
                    event_observed_A=high["OS"],
                    event_observed_B=low["OS"]
                )
                cph = CoxPHFitter()
                cph.fit(
                    sub[["OS.time", "OS", "High_Group"]],
                    duration_col="OS.time",
                    event_col="OS"
                )
                rows.append({
                    "Cancer": cancer,
                    "Patients": len(sub),
                    "Deaths": int(sub["OS"].sum()),
                    "High deaths": int(high["OS"].sum()),
                    "Low deaths": int(low["OS"].sum()),
                    "Log-rank P": float(lr.p_value),
                    "HR": float(cph.hazard_ratios_["High_Group"]),
                    "HR lower": float(np.exp(
                        cph.confidence_intervals_.loc["High_Group"].iloc[0]
                    )),
                    "HR upper": float(np.exp(
                        cph.confidence_intervals_.loc["High_Group"].iloc[1]
                    )),
                    "Cox P": float(cph.summary.loc["High_Group", "p"])
                })
            except Exception as exc:
                errors.append({
                    "Cancer": cancer,
                    "Reason": f"模型估計失敗：{str(exc)[:250]}"
                })

        if rows:
            result = pd.DataFrame(rows)
            result["Log-rank FDR"] = benjamini_hochberg(
                result["Log-rank P"]
            )
            result["Cox FDR"] = benjamini_hochberg(result["Cox P"])
            result = result.sort_values("Cox FDR")
            st.session_state["role_survival_results"] = result
        else:
            st.session_state["role_survival_results"] = None

        st.session_state["role_survival_errors"] = pd.DataFrame(errors)

    result = st.session_state.get("role_survival_results")
    errors = st.session_state.get("role_survival_errors")

    if result is not None:
        st.subheader("Pan-Cancer survival summary")
        st.dataframe(
            result.style.format({
                "Log-rank P": "{:.4g}",
                "HR": "{:.3f}",
                "HR lower": "{:.3f}",
                "HR upper": "{:.3f}",
                "Cox P": "{:.4g}",
                "Log-rank FDR": "{:.4g}",
                "Cox FDR": "{:.4g}"
            }),
            use_container_width=True,
            hide_index=True
        )

        significant = result[result["Cox FDR"] < 0.05]
        st.metric("Cancer types with Cox FDR < 0.05", len(significant))
        if not significant.empty:
            st.write("Cox FDR < 0.05 的癌症")
            st.dataframe(
                significant, use_container_width=True, hide_index=True
            )

        st.download_button(
            "下載 Role Assessment 生存證據 CSV",
            result.to_csv(index=False).encode("utf-8-sig"),
            file_name="TARBP2_Role_Assessment_Survival.csv",
            mime="text/csv"
        )
    elif "role_survival_results" in st.session_state:
        st.warning("沒有癌症種類成功完成生存模型，請查看下方失敗原因。")

    if errors is not None and not errors.empty:
        with st.expander("查看未納入分析的癌症及原因"):
            st.dataframe(errors, use_container_width=True, hide_index=True)
except Exception as exc:
    st.warning(f"生存證據目前無法載入：{exc}")


# ---------------------------------------------------------
# 4. RBP / STRING network evidence
# ---------------------------------------------------------
st.divider()
st.header("4. RBP / STRING Network Evidence")

st.markdown("""
STRING 的 combined score 是功能關聯的信心分數，不是結合親和力、
相關係數或因果機率。Functional association 也不等於已證實的直接物理互作。
""")

network_path = find_existing([
    "data/STRING/TARBP2_STRING_network.csv",
    "data/STRING/TARBP2_network.csv",
    "data/RBP/TARBP2_STRING_network.csv",
    "data/TARBP2_STRING_network.csv",
    "data/GEO/TARBP2_STRING_network.csv"
])
network = read_csv_safe(network_path)

uploaded_network = st.file_uploader(
    "若網絡 CSV 未自動找到，可在此上傳 STRING network CSV",
    type=["csv"],
    key="role_network_upload"
)
if uploaded_network is not None:
    try:
        network = pd.read_csv(uploaded_network)
        st.success("已載入上傳的 STRING network CSV。")
    except Exception as exc:
        st.error(f"讀取上傳檔案失敗：{exc}")

if network is None:
    st.info(
        "尚未找到可供本頁讀取的 STRING network CSV。"
        "你仍可在 RBP Network 頁面查看網絡；若要在本頁整合，"
        "請下載該頁的 CSV 後上傳。"
    )
else:
    st.write(f"載入網絡資料：**{len(network)}** 筆連線")
    st.dataframe(network.head(30), use_container_width=True, hide_index=True)
    score_col = next(
        (
            c for c in network.columns
            if c.lower() in ["combined_score", "combined score", "score"]
        ),
        None
    )
    if score_col:
        scores = pd.to_numeric(
            network[score_col], errors="coerce"
        ).dropna()
        if not scores.empty:
            st.metric("Highest STRING score", f"{scores.max():.3f}")
            st.caption(
                "請確認匯出分數尺度：STRING API 常見 combined score 為 0–1000；"
                "若已正規化則可能為 0–1。此數值僅代表資料中的最高分。"
            )
    st.download_button(
        "下載目前載入的網絡資料",
        network.to_csv(index=False).encode("utf-8-sig"),
        file_name="TARBP2_RBP_Network_Evidence.csv",
        mime="text/csv"
    )


# ---------------------------------------------------------
# 5. Integrated interpretation
# ---------------------------------------------------------
st.divider()
st.header("5. Integrated Evidence Interpretation")

st.markdown("""
| 證據類型 | 可以支持的判讀 | 不能直接推論的事項 |
|---|---|---|
| Pan-cancer expression | TARBP2 在不同癌症樣本中的表現分布 | 表現較高就代表促癌 |
| Tumor vs Normal | 各癌症腫瘤與正常樣本的表現差異 | 差異由 TARBP2 導致 |
| Overall survival | TARBP2 表現與 OS 的統計關聯 | TARBP2 是因果因素或已證實的獨立預後因子 |
| STRING network | TARBP2 與候選 RBP/RNAi 蛋白的功能關聯線索 | 所有連線均為直接物理互作 |
""")

st.subheader("目前可採用的研究結論")
st.success("""
本研究利用 TCGA 公開資料分析 TARBP2 在不同癌症中的表現分布、
腫瘤與正常組織差異，以及其與整體存活率的統計關聯，並結合
STRING 蛋白質功能關聯網絡，探討 TARBP2 在癌症相關 RNA 調控中的潛在角色。
上述結果屬於觀察性資料所提供的關聯證據，可作為後續研究假說，
尚不足以單獨證明 TARBP2 對腫瘤具有促進或抑制的因果作用。
""")

st.warning("""
研究限制：
1. TCGA 為觀察性資料，無法單獨建立因果關係。
2. Pan-Cancer Cox 為單變項分析，未控制年齡、性別、分期等潛在混雜因子。
3. Tumor vs Normal 依癌症種類分開比較；配對樣本使用配對檢定，仍未校正臨床共變項。
4. 病患配對依 TCGA sample barcode 前 12 個字元辨識，須確認原始資料 barcode 格式正確。
5. STRING 的功能關聯分數不等於直接結合或實驗驗證。
6. 各癌症樣本數、事件數及資料完整度不同，結果需依各癌症的分析樣本數解讀。
""")

st.caption(
    "TARBP2 Role Assessment | Research and educational use only. "
    "Not intended for clinical decision-making."
)
