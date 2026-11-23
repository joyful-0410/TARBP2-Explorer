import streamlit as st

st.set_page_config(
    page_title="TRBP2 & RNAi",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 TRBP2 & RNAi")
st.markdown(
    """
    本頁介紹 **TARBP2 基因所編碼的 TRBP2 蛋白**，以及它在小 RNA
    生合成與 RNA interference（RNAi）途徑中的角色。
    """
)

st.info(
    "名詞提醒：TARBP2 是基因名稱；TRBP2 是其蛋白質產物的常用名稱。"
    "部分文獻也會使用 TRBP（TAR RNA-binding protein）等稱呼，閱讀文獻時需留意上下文。"
)

tab1, tab2, tab3, tab4 = st.tabs(
    ["RNAi 基本概念", "miRNA 生合成路徑", "TRBP2 的角色", "重要名詞與文獻"]
)

with tab1:
    st.header("1. 什麼是 RNA interference（RNAi）？")
    st.write(
        """
        RNA interference（RNAi，RNA 干擾）是一種由小型 RNA 分子引導的基因表現調控機制。
        小 RNA 與 Argonaute（AGO）蛋白等組成 RNA-induced silencing complex（RISC），
        再由其中的 guide strand 引導複合體辨識目標 mRNA，降低目標基因的表現。
        """
    )

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("miRNA")
        st.write(
            """
            - 通常由細胞自身的 miRNA 基因轉錄產生。
            - 常與目標 mRNA 部分互補。
            - 常造成翻譯抑制、mRNA 去腺苷化及後續降解。
            """
        )
    with col2:
        st.subheader("siRNA")
        st.write(
            """
            - 可由外源雙股 RNA 或人工設計的 siRNA 產生。
            - 與目標 mRNA 通常具有高度互補性。
            - 在典型情況下可引導 AGO2 切割目標 mRNA。
            """
        )

    st.caption(
        "miRNA 與 siRNA 都可利用 Argonaute/RISC 進行基因沉默；"
        "兩者的來源、序列配對特性及典型作用方式有所不同。"
    )

with tab2:
    st.header("2. miRNA 的典型生合成路徑")
    st.markdown(
        """
        **典型 miRNA 路徑可概略分成以下步驟：**

        1. **轉錄**：miRNA 基因通常由 RNA polymerase II 轉錄形成 pri-miRNA。
        2. **核內加工**：Drosha 與 DGCR8 組成的 Microprocessor complex 將 pri-miRNA
           加工成較短的 hairpin precursor（pre-miRNA）。
        3. **輸出至細胞質**：Exportin-5/Ran-GTP 協助 pre-miRNA 從細胞核輸出。
        4. **Dicer 加工**：Dicer 將 pre-miRNA 切割成短的 miRNA duplex；
           TRBP2 可與 Dicer 形成複合體並參與此步驟。
        5. **RISC loading**：miRNA duplex 被載入含 Argonaute 的沉默複合體，
           通常保留其中一股作為 guide strand。
        6. **目標辨識與沉默**：guide strand 引導複合體辨識目標 mRNA，
           造成翻譯抑制或 mRNA 不穩定及降解。
        """
    )

    st.warning(
        "這裡呈現的是典型 miRNA 生合成流程。不同 miRNA 可能存在非典型生合成路徑，"
        "因此不應將所有小 RNA 都視為完全相同的加工流程。"
    )

with tab3:
    st.header("3. TRBP2 在 RNAi 路徑中的角色")
    st.write(
        """
        TRBP2（TAR RNA-binding protein 2）是由 **TARBP2** 基因編碼的雙股 RNA
        結合蛋白。在典型 miRNA 加工過程中，TRBP2 可與 Dicer 互相作用形成複合體，
        參與 pre-miRNA 的辨識與加工，並可能協助 Dicer 與 Argonaute 之間的小 RNA
        傳遞及沉默複合體組裝。不同研究系統對 TRBP2 的必要性與作用程度並非完全一致。
        """
    )

    st.subheader("TRBP2 與 Dicer")
    st.markdown(
        """
        - **Dicer**：RNase III 家族核酸內切酶，可將 pre-miRNA 等雙股 RNA 前驅物
          加工成短小 RNA duplex。
        - **TRBP2**：具有雙股 RNA 結合能力，可與 Dicer 形成複合體，
          協助小 RNA 加工過程。
        - **共同作用**：Dicer–TRBP2 複合體可促進部分 miRNA 前驅物的辨識與加工，並與後續沉默複合體形成相關。
        """
    )

    st.subheader("簡化關係圖")
    st.code(
        "pre-miRNA\n"
        "    ↓\n"
        "Dicer  +  TRBP2\n"
        "    ↓\n"
        "miRNA duplex\n"
        "    ↓\n"
        "Argonaute / RISC loading\n"
        "    ↓\n"
        "guide strand 辨識目標 mRNA\n"
        "    ↓\n"
        "基因表現受到抑制",
        language="text"
    )

    st.subheader("與本專題 TARBP2 的連結")
    st.write(
        """
        本專題研究的 TARBP2 是 TRBP2 蛋白的編碼基因。因此，分析 TARBP2 的表現量時，
        可以從「小 RNA 加工與 RNA 調控」的生物學背景切入，討論它可能與細胞基因表現
        調控之間的關聯。
        """
    )
    st.warning(
        "TARBP2 表現量與癌症臨床結果之間若出現統計關聯，不能單憑相關性就推論"
        "RNAi 活性一定改變，也不能直接推論 TRBP2 透過某一特定 miRNA 造成癌症表型。"
        "這些機制仍需要額外的分子實驗或文獻證據支持。"
    )

    st.info(
        "研究證據補充：早期研究支持 TRBP 可協助 Dicer 加工小 RNA，並與 Ago2/RISC 組裝相關；"
        "後續 TRBP knockout 研究則指出，TRBP 對整體 miRNA 豐度或 Argonaute loading 的影響"
        "未必在所有細胞條件下都能觀察到。因此較嚴謹的說法是「TRBP2 為 Dicer 的調控因子／"
        "輔因子之一」，不宜寫成所有 miRNA 加工都絕對依賴 TRBP2。"
    )

with tab4:
    st.header("4. 重要名詞整理")
    terms = [
        ("TARBP2", "編碼 TRBP2 蛋白的基因名稱。"),
        ("TRBP2", "TARBP2 基因所編碼的雙股 RNA 結合蛋白。"),
        ("RNAi", "由小 RNA 引導、降低目標基因表現的調控機制。"),
        ("miRNA", "通常由細胞內源性產生的小型非編碼 RNA，可調控目標 mRNA。"),
        ("siRNA", "小干擾 RNA，常與目標 mRNA 高度互補並引導其沉默。"),
        ("Dicer", "RNase III 家族核酸內切酶，參與小 RNA 前驅物加工。"),
        ("Argonaute（AGO）", "RISC 的核心蛋白家族，可結合 guide RNA 並參與目標辨識。"),
        ("RISC", "RNA-induced silencing complex，由 guide RNA 與蛋白質組成，負責目標辨識及沉默。"),
        ("Gene silencing", "基因沉默，指基因表現受到抑制；RNAi 是其中一種機制。"),
        ("pre-miRNA", "由 pri-miRNA 加工而成、具有髮夾結構的 miRNA 前驅物。"),
        ("miRNA duplex", "Dicer 加工後形成的短雙股 RNA，後續可進入 RISC loading。"),
    ]

    for term, definition in terms:
        with st.expander(term):
            st.write(definition)

    st.subheader("閱讀文獻時要記得")
    st.markdown(
        """
        - **基因與蛋白質不要混用**：TARBP2 指基因，TRBP2 指蛋白質。
        - **RNAi 不等於只有 miRNA**：siRNA 也是常見的小 RNA 沉默途徑。
        - **TRBP2 不等於 Dicer**：兩者是不同蛋白，但可形成複合體並共同參與小 RNA 加工。
        - **路徑示意不等於已證實本專題機制**：若要主張 TARBP2 在特定癌症中透過某條 RNAi 路徑發揮作用，需要相應證據。
        """
    )

    st.subheader("核心參考文獻")
    st.markdown(
        """
        1. Chendrimada et al. (2005). **TRBP recruits the Dicer complex to Ago2 for microRNA processing and gene silencing.** *Nature*, 436, 740–744.  
           [PubMed](https://pubmed.ncbi.nlm.nih.gov/15973356/) · [DOI: 10.1038/nature03868](https://doi.org/10.1038/nature03868)

        2. Wilson et al. (2010). **Characterization of the TRBP domain required for Dicer interaction and function in RNA interference.** *BMC Biochemistry*, 11, 38.  
           [PubMed](https://pubmed.ncbi.nlm.nih.gov/19422693/)

        3. Kim et al. (2014). **Deletion of human tarbp2 reveals cellular microRNA targets and cell-cycle function of TRBP.** *Cell Reports*.  
           [PubMed](https://pubmed.ncbi.nlm.nih.gov/25437560/)

        4. Fareh et al. (2016). **TRBP ensures efficient Dicer processing of precursor microRNA in RNA-crowded environments.** *Nature Communications*, 7, 13694.  
           [文章全文](https://www.nature.com/articles/ncomms13694)
        """
    )
    st.caption(
        "閱讀建議：第 1、2 篇可了解 Dicer–TRBP2 互作與 RNAi 功能；"
        "第 3 篇提醒 TRBP2 在 miRNA 生合成中的作用具有情境依賴性；"
        "第 4 篇補充 TRBP 如何協助 Dicer 選擇與加工 pre-miRNA。"
    )

st.divider()
st.caption(
    "本頁為 RNAi 與 TRBP2 的基礎知識介紹。內容用於研究背景與教學理解，"
    "不代表本專題已直接驗證上述所有分子機制。"
)
