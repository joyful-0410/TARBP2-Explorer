from pathlib import Path

import streamlit as st

st.set_page_config(
    page_title="TARBP2 Protein Structure",
    page_icon="🧬",
    layout="wide",
)

# Optional 3D-viewer dependencies: install with
# pip install py3Dmol stmol
try:
    import py3Dmol
    from stmol import showmol
    VIEWER_AVAILABLE = True
except ImportError:
    VIEWER_AVAILABLE = False


# This page is stored in the project's pages/ directory.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
PDB_PATH = PROJECT_ROOT / "data" / "structure" / "TARBP2.pdb"
IMAGE_PATH = PROJECT_ROOT / "images" / "trbp2_structure.png"


# Literature-based residue ranges for human TRBP2 (366 aa).
# dsRBD1: 31–96; dsRBD2: 160–226; C4: 298–366.
# The broader Medipal region is approximately 237–366.
DOMAINS = {
    "Whole Protein": None,
    "dsRBD1": (31, 96),
    "dsRBD2": (160, 226),
    "C4 / Dicer-binding domain": (298, 366),
}

st.title("🧬 TARBP2 / TRBP2 Protein Structure")
st.write(
    """
    本頁介紹人類 TARBP2 基因所編碼的 TRBP2 蛋白質（366 aa），
    以 AlphaFold 預測結構作為視覺化參考，並整理其雙股 RNA 結合區域
    與 C-terminal protein-interaction domain 的功能。
    """
)

st.info(
    """
    **結構解讀提醒：** AlphaFold 提供的是蛋白質結構預測，不是本研究自行
    取得的實驗結構。圖中的區域標示依據已發表文獻的胺基酸位置；
    若 PDB 檔案採用不同的殘基編號或包含額外標籤，顏色標示可能需要再校正。
    """
)

st.divider()
st.header("1. TRBP2 的主要結構區域")

c1, c2, c3 = st.columns(3)

with c1:
    st.markdown("### dsRBD1")
    st.success(
        """
        **Residues 31–96**

        第一個 double-stranded RNA-binding domain，
        參與雙股 RNA 的辨識與結合。
        """
    )

with c2:
    st.markdown("### dsRBD2")
    st.warning(
        """
        **Residues 160–226**

        第二個 dsRNA-binding domain。
        與 dsRBD1 協同結合 RNA；文獻指出其 dsRNA-binding
        活性通常較強。
        """
    )

with c3:
    st.markdown("### C4 domain")
    st.info(
        """
        **Residues 298–366**

        C-terminal protein-interaction domain，
        參與 TRBP2 與 Dicer 等蛋白質的交互作用；
        此區域本身不應描述成第三個典型 RNA-binding domain。
        """
    )

st.markdown(
    """
    **補充：** Medipal region 約為 residues 237–366，範圍比 C4 domain
    更廣；文獻指出 C4 domain 足以支持 TRBP–Dicer interaction。
    """
)

st.divider()
st.header("2. AlphaFold 預測結構")

if IMAGE_PATH.exists():
    st.image(
        str(IMAGE_PATH),
        caption="AlphaFold predicted structure of human TRBP2 / TARBP2",
        use_container_width=True,
    )
else:
    st.warning(
        f"找不到結構示意圖片：{IMAGE_PATH.relative_to(PROJECT_ROOT)}。"
        "若有需要，請確認 images 資料夾內是否存在 trbp2_structure.png。"
    )

st.caption(
    "AlphaFold predicted structure；顏色示意應與下方互動式檢視器的區域定義一致。"
)

st.divider()
st.header("3. Interactive 3D Viewer")

domain = st.selectbox("選擇要檢視的區域", list(DOMAINS.keys()))

if not PDB_PATH.exists():
    st.error(
        f"找不到 PDB 檔案：{PDB_PATH.relative_to(PROJECT_ROOT)}。"
        "請確認 data/structure/TARBP2.pdb 已放在專題資料夾中。"
    )
elif not VIEWER_AVAILABLE:
    st.warning(
        "尚未安裝 3D 檢視器套件。請在終端機執行："
        "`pip install py3Dmol stmol`，再重新啟動 Streamlit。"
    )
else:
    try:
        pdb_data = PDB_PATH.read_text(encoding="utf-8", errors="replace")

        viewer = py3Dmol.view(width=900, height=620)
        viewer.addModel(pdb_data, "pdb")
        viewer.setStyle({"cartoon": {"color": "lightgray"}})

        selected_range = DOMAINS[domain]

        if selected_range is None:
            viewer.setStyle({"cartoon": {"color": "spectrum"}})
        else:
            start_residue, end_residue = selected_range
            residue_numbers = list(range(start_residue, end_residue + 1))
            highlight_color = {
                "dsRBD1": "green",
                "dsRBD2": "orange",
                "C4 / Dicer-binding domain": "purple",
            }[domain]

            viewer.setStyle(
                {"resi": residue_numbers},
                {"cartoon": {"color": highlight_color}},
            )

        viewer.zoomTo()
        showmol(viewer, height=620, width=900)

    except Exception as exc:
        st.error(f"讀取或呈現 PDB 結構時發生錯誤：{exc}")

st.divider()
st.header("4. 區域功能說明")

if domain == "Whole Protein":
    st.subheader("Whole Protein")
    st.info(
        """
        **Whole Protein | Residues 1–366**

        TRBP2 全長約 366 個胺基酸，包含兩個典型 dsRNA-binding domains，
        以及 C-terminal C4 domain。兩個 dsRBD 參與 RNA 結合；
        C4 domain 則主要參與蛋白質交互作用，包括與 Dicer 的結合。
        """
    )

elif domain == "dsRBD1":
    st.subheader("dsRBD1 | Residues 31–96")
    st.success(
        """
        **Residues 31–96**

        - 屬於典型 double-stranded RNA-binding domain。
        - 參與雙股 RNA 的辨識與結合。
        - 與 dsRBD2 協同作用，可提高 TRBP2 對 RNA 的結合能力。
        """
    )

elif domain == "dsRBD2":
    st.subheader("dsRBD2 | Residues 160–226")
    st.warning(
        """
        **Residues 160–226**

        - 第二個典型 dsRNA-binding domain。
        - 具有與 dsRBD1 相似的 RNA-binding fold。
        - 文獻指出 dsRBD2 的 dsRNA-binding 活性通常較強，
          並與其 KR-helix motif 等結構特徵有關。
        """
    )

else:
    st.subheader("C4 domain | Residues 298–366")
    st.info(
        """
        **Residues 298–366**

        - 位於 TRBP2 的 C-terminal。
        - 是 TRBP2 與 Dicer 交互作用的重要區域。
        - 雖具有部分 dsRBD-like structural features，
          但不應將 C4 描述為典型的 RNA-binding domain。
        - 更廣的 Medipal region 約為 residues 237–366。
        """
    )

st.divider()
st.header("5. Biological significance")
st.write(
    """
    TRBP2 可透過兩個 dsRBD 與雙股 RNA 結合，並透過 C-terminal C4
    domain 與 Dicer 互動。這些結構特徵使 TRBP2 能參與小 RNA 加工及
    RNA interference 相關路徑。此頁呈現的是已發表文獻與預測結構的
    整合解讀，不能單憑結構圖推論 TARBP2 在特定癌症中的因果作用。
    """
)

st.header("6. References")
st.markdown(
    """
    1. **The Multiple Functions of TRBP, at the Hub of Cell Responses to Viruses, Stress, and Cancer.**
       *Microbiology and Molecular Biology Reviews* (2012).
    2. **Characterization of the TRBP domain required for dicer interaction and function in RNA interference.**
       *Journal of Biological Chemistry* (2009).
    3. **Structures of the first and second double-stranded RNA-binding domains of human TAR RNA-binding protein.**
       *Journal of Molecular Biology* (2010).
    """
)
