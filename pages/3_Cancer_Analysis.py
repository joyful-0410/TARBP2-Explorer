import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================
# TARBP2 Pan-Cancer Expression Analysis
# ==========================================

st.set_page_config(
    page_title="TARBP2 Cancer Analysis",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 TARBP2 Pan-Cancer Expression Analysis")

st.markdown("""
本模組利用 TCGA Pan-Cancer RNA-seq 資料，
分析 TARBP2 在不同癌症中的基因表現量。
""")

st.divider()

# ==========================================
# 1. Load data
# ==========================================

data_file = "data/TCGA/TARBP2_primary_tumor.csv"
summary_file = "data/TCGA/TARBP2_cancer_summary.csv"

try:
    df = pd.read_csv(data_file)
    summary = pd.read_csv(summary_file)

except FileNotFoundError:
    st.error("找不到 TCGA 分析資料，請先執行 merge_tcga.py")
    st.stop()

# ==========================================
# 2. Data cleaning
# ==========================================

df["TARBP2_expression"] = pd.to_numeric(
    df["TARBP2_expression"],
    errors="coerce"
)

df = df.dropna(
    subset=["TARBP2_expression", "_primary_disease"]
)

# ==========================================
# 3. Overview
# ==========================================

col1, col2, col3 = st.columns(3)

col1.metric(
    "Total Tumor Samples",
    f"{len(df):,}"
)

col2.metric(
    "Cancer Types",
    f"{df['_primary_disease'].nunique()}"
)

col3.metric(
    "Gene",
    "TARBP2"
)

st.divider()

# ==========================================
# 4. Cancer summary
# ==========================================

summary = summary.sort_values(
    "median",
    ascending=False
)

st.subheader("1. Pan-Cancer Expression Ranking")

st.dataframe(
    summary.rename(columns={
        "_primary_disease": "Cancer Type",
        "sample_count": "Sample Count",
        "mean": "Mean Expression",
        "median": "Median Expression",
        "std": "Standard Deviation"
    }),
    use_container_width=True,
    hide_index=True
)

# ==========================================
# 5. Box Plot
# ==========================================

st.divider()

st.subheader("2. TARBP2 Expression Distribution")

st.write(
    "比較不同癌症中的 TARBP2 RNA expression distribution。"
)

cancer_order = (
    df.groupby("_primary_disease")["TARBP2_expression"]
    .median()
    .sort_values(ascending=False)
    .index.tolist()
)

plot_data = [
    df.loc[
        df["_primary_disease"] == cancer,
        "TARBP2_expression"
    ].values
    for cancer in cancer_order
]

fig, ax = plt.subplots(figsize=(12, 14))

ax.boxplot(
    plot_data,
    vert=False,
    showfliers=False
)

ax.set_yticks(range(1, len(cancer_order) + 1))
ax.set_yticklabels(cancer_order)

ax.set_xlabel("TARBP2 Expression (log2 normalized value + 1)")
ax.set_ylabel("Cancer Type")

ax.set_title(
    "TARBP2 Pan-Cancer Expression Distribution"
)

ax.invert_yaxis()

plt.tight_layout()

st.pyplot(fig)

plt.close(fig)

# ==========================================
# 6. Cancer selection
# ==========================================

st.divider()

st.subheader("3. Individual Cancer Analysis")

selected_cancer = st.selectbox(
    "Select Cancer Type",
    cancer_order
)

selected_data = df[
    df["_primary_disease"] == selected_cancer
]

col1, col2, col3 = st.columns(3)

col1.metric(
    "Sample Count",
    len(selected_data)
)

col2.metric(
    "Mean Expression",
    f"{selected_data['TARBP2_expression'].mean():.3f}"
)

col3.metric(
    "Median Expression",
    f"{selected_data['TARBP2_expression'].median():.3f}"
)

fig2, ax2 = plt.subplots(figsize=(8, 4))

ax2.hist(
    selected_data["TARBP2_expression"],
    bins=25
)

ax2.set_title(
    f"TARBP2 Expression in {selected_cancer}"
)

ax2.set_xlabel("Expression")
ax2.set_ylabel("Frequency")

st.pyplot(fig2)

plt.close(fig2)

# ==========================================
# 7. Download data
# ==========================================

st.divider()

st.subheader("4. Download Analysis Data")

st.download_button(
    label="Download TARBP2 Expression CSV",
    data=df.to_csv(index=False).encode("utf-8-sig"),
    file_name="TARBP2_PanCancer_Expression.csv",
    mime="text/csv"
)

st.success("TARBP2 Pan-Cancer Analysis loaded successfully.")

st.caption(
    "Data source: UCSC Xena TCGA Pan-Cancer Atlas. "
    "Descriptive analysis of primary tumor samples. "
    "Expression differences do not establish causality."
)