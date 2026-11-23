import streamlit as st

st.set_page_config(
    page_title="TARBP2 Cancer Bioinformatics Platform",
    page_icon="🧬",
    layout="wide",
)

# =========================================================
# TARBP2 Cancer Bioinformatics Platform — Home
# 本檔只負責首頁；分析功能由 pages/ 資料夾中的頁面提供。
# Streamlit 會依 pages/ 中的 .py 檔建立多頁導覽。
# =========================================================

st.title("🧬 TARBP2 Cancer Bioinformatics Platform")
st.subheader(
    "Public Cancer Data × Survival Analysis × RBP Network × Protein Structure"
)

st.markdown(
    """
    本平台以 **TARBP2 / TRBP2** 為研究核心，整合公開癌症資料、
    腫瘤與正常組織表現比較、生存分析、RNA-binding protein（RBP）/
    RNA interference（RNAi）功能網絡，以及蛋白質結構資訊，
    探索 TARBP2 在不同癌症情境中的表現特徵與可能的生物學意義。
    """
)

st.info(
    """
    **研究定位與範圍**

    本專題以公開資料庫與生物資訊工具進行計算分析，屬於探索性、
    觀察性研究。研究內容不包含新進行的濕實驗驗證；
    平台結果用於提出及整理研究假說，不直接代表因果關係，
    也不等同於臨床診斷或治療建議。
    """
)

st.divider()

st.header("研究流程")
st.markdown(
    """
    **公開資料整理 → TARBP2 表現分析 → 腫瘤／正常比較
    → 生存與預後分析 → RBP / RNAi 網絡探索
    → 蛋白質結構資訊整合 → 綜合角色評估**
    """
)

st.divider()

st.header("平台主要模組")

col1, col2 = st.columns(2)

with col1:
    st.markdown("### 🧬 TARBP2 & RNAi")
    st.write(
        "介紹 TARBP2 / TRBP2、Dicer、miRNA maturation 與 "
        "RNA interference（RNAi）機制，作為理解研究背景的文獻基礎。"
    )

    st.markdown("### 📊 Cancer Analysis")
    st.write(
        "探索 TARBP2 在不同癌症中的表現分布，整理公開癌症資料庫的分析結果。"
    )

    st.markdown("### 🔗 RBP Network")
    st.write(
        "使用 STRING 等公開資源探索 TARBP2 與 DICER1、AGO2、DROSHA、"
        "DGCR8、IGF2BP1、ELAVL1 等候選蛋白質的功能關聯。"
    )

with col2:
    st.markdown("### 📈 Survival Analysis")
    st.write(
        "以 Kaplan–Meier、Log-rank 與 Cox regression 評估 TARBP2 表現"
        "與 overall survival（OS）的統計關聯。"
    )

    st.markdown("### 🧩 Protein Structure")
    st.write(
        "整理 TARBP2 編碼之 TRBP2 蛋白質結構資訊，介紹 dsRNA-binding "
        "domains（dsRBDs）及 C-terminal 區域的可能功能。"
    )

    st.markdown("### 🧭 Role Assessment")
    st.write(
        "整合表現量、腫瘤／正常比較、生存分析與 RBP network 證據，"
        "進行 TARBP2 候選角色的探索性評估。"
    )

st.divider()

st.header("研究問題")
st.write(
    """
    TARBP2 基因編碼 TRBP2 蛋白質。TRBP2 具有雙股 RNA 結合區域，
    並與 Dicer 及小 RNA 加工路徑相關。本平台藉由公開資料分析與
    已發表的分子功能資訊，探索 TARBP2 在不同癌症中的表現、
    預後關聯，以及其與 RNA-binding proteins 的功能網絡線索。
    """
)

st.header("資料解讀與研究限制")
st.warning(
    """
    - 公開資料多屬觀察性資料，統計關聯不能直接解讀為因果關係。
    - 腫瘤與正常樣本的比較需考慮組織來源、癌症種類及樣本是否配對。
    - 生存分析可能受到年齡、性別、分期等臨床變項及模型假設影響。
    - 跨癌症比較需留意樣本數、臨床組成及多重比較問題。
    - STRING functional association 代表功能關聯線索，不等於已證實直接物理交互作用。
    - AlphaFold 結構屬預測模型，不等同於本研究自行取得的實驗結構。
    - 本平台不包含濕實驗驗證，結論應限於公開資料與計算分析所支持的範圍。
    - 平台結果屬研究探索用途，不可作為個人診斷或治療建議。
    """
)

st.caption(
    "TARBP2 Cancer Bioinformatics Platform | Public-data-based computational research"
)
