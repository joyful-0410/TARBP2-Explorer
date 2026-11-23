import io
import requests
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

try:
    import networkx as nx
except ImportError:
    nx = None


st.set_page_config(
    page_title="TARBP2 RBP Network",
    page_icon="🔗",
    layout="wide"
)

st.title("🔗 TARBP2-associated RBP / RNAi Network")
st.markdown(
    """
    本頁使用 **STRING database** 查詢人類 TARBP2 與候選 RNA-binding /
    RNAi-related proteins 之間的蛋白質功能關聯網絡，並呈現互作分數與證據通道。
    """
)

st.info(
    "解讀提醒：STRING 的 functional network 代表蛋白質之間已知或預測的功能關聯，"
    "不等於每一對蛋白質都已證實直接物理結合。若選擇 physical network，"
    "才會聚焦於物理互作證據；即使如此，仍須依資料庫證據解讀。"
)

STRING_API = "https://version-12-5.string-db.org/api"
CALLER_IDENTITY = "TARBP2_Explorer"

# Candidate panel: TARBP2 is the query protein; the remaining proteins are
# biological-context candidates, not assumed interactors.
DEFAULT_PROTEINS = [
    "TARBP2",
    "DICER1",
    "AGO2",
    "DROSHA",
    "DGCR8",
    "IGF2BP1",
    "ELAVL1",
]

PROTEIN_NOTES = {
    "TARBP2": "本專題核心基因；編碼 TRBP2 雙股 RNA 結合蛋白。",
    "DICER1": "RNase III family endoribonuclease；參與小 RNA 前驅物加工。",
    "AGO2": "Argonaute family protein；參與 RISC 與小 RNA 導引的基因沉默。",
    "DROSHA": "Microprocessor complex 的 RNase III enzyme；參與 pri-miRNA 核內加工。",
    "DGCR8": "Microprocessor complex component；協助辨識 pri-miRNA。",
    "IGF2BP1": "RNA-binding protein；參與 RNA 定位、穩定性及翻譯調控。",
    "ELAVL1": "RNA-binding protein；可影響多種 mRNA 的穩定性與命運。",
}


@st.cache_data(ttl=86400, show_spinner=False)
def map_proteins(protein_names):
    url = f"{STRING_API}/tsv/get_string_ids"
    response = requests.get(
        url,
        params={
            "identifiers": "\r".join(protein_names),
            "species": 9606,
            "caller_identity": CALLER_IDENTITY,
        },
        timeout=30,
    )
    response.raise_for_status()
    if not response.text.strip():
        raise ValueError("STRING 沒有回傳蛋白質 ID 對應結果。")

    mapped = pd.read_csv(io.StringIO(response.text), sep="\t")
    # STRING API v12 returns queryIndex (not queryItem) to identify the
    # position of each submitted identifier. Reconstruct the original query
    # name from that index so the rest of the app uses a stable queryItem field.
    required = {"queryIndex", "stringId", "preferredName"}
    if not required.issubset(mapped.columns):
        raise ValueError(
            "STRING ID 對應結果缺少必要欄位。"
            f"需要欄位：{sorted(required)}；實際欄位：{list(mapped.columns)}"
        )

    mapped["queryIndex"] = pd.to_numeric(mapped["queryIndex"], errors="coerce")
    mapped = mapped.dropna(subset=["queryIndex", "stringId", "preferredName"]).copy()
    mapped["queryIndex"] = mapped["queryIndex"].astype(int)

    # STRING queryIndex is zero-based and follows the submitted identifier order.
    mapped["queryItem"] = mapped["queryIndex"].apply(
        lambda idx: protein_names[idx] if 0 <= idx < len(protein_names) else None
    )
    mapped = mapped.dropna(subset=["queryItem"])
    return mapped


@st.cache_data(ttl=86400, show_spinner=False)
def get_network(string_ids, required_score, network_type):
    url = f"{STRING_API}/tsv/network"
    response = requests.get(
        url,
        params={
            "identifiers": "\r".join(string_ids),
            "species": 9606,
            "required_score": int(required_score),
            "network_type": network_type,
            "add_nodes": 0,
            "caller_identity": CALLER_IDENTITY,
        },
        timeout=45,
    )

    if response.status_code == 404:
        return pd.DataFrame()
    response.raise_for_status()

    if not response.text.strip():
        return pd.DataFrame()

    return pd.read_csv(io.StringIO(response.text), sep="\t")


with st.sidebar:
    st.header("Network settings")
    network_type_label = st.selectbox(
        "網絡類型",
        ["Functional association", "Physical interaction"],
        help=(
            "Functional：功能關聯網絡；Physical：物理互作網絡。"
            "兩者的生物學意義不同。"
        ),
    )
    network_type = (
        "functional"
        if network_type_label == "Functional association"
        else "physical"
    )

    score_threshold = st.slider(
        "最低 STRING combined score",
        min_value=150,
        max_value=900,
        value=400,
        step=50,
        help="STRING combined score 範圍為 0–1000；提高門檻會留下較高分數的連線。",
    )

    selected_proteins = st.multiselect(
        "候選蛋白質",
        options=DEFAULT_PROTEINS,
        default=DEFAULT_PROTEINS,
        help="這是本專題設定的候選清單；不代表這些蛋白質彼此已確定互作。",
    )

    run_query = st.button("查詢 STRING 網絡", type="primary", use_container_width=True)


st.subheader("1. 本次查詢設定")
st.write(
    f"物種：**Homo sapiens（9606）**　|　"
    f"網絡：**{network_type_label}**　|　"
    f"最低分數：**{score_threshold}/1000**"
)
st.caption(
    "本頁只查詢所選候選蛋白之間的網絡（add_nodes=0），"
    "不會自動加入未選取的鄰近蛋白。"
)

if not selected_proteins:
    st.warning("請至少選擇一個候選蛋白質。")
    st.stop()

if run_query or "rbp_network_result" not in st.session_state:
    with st.spinner("正在連線 STRING、對應蛋白質 ID 並取得網絡資料……"):
        try:
            mapping = map_proteins(selected_proteins)
            mapped_names = set(mapping["queryItem"].astype(str))
            missing = [p for p in selected_proteins if p not in mapped_names]

            if mapping.empty:
                st.error("STRING 無法對應所選蛋白質。")
                st.stop()

            ids = mapping["stringId"].dropna().astype(str).drop_duplicates().tolist()
            network = get_network(ids, score_threshold, network_type)

            st.session_state["rbp_network_result"] = {
                "mapping": mapping,
                "network": network,
                "missing": missing,
                "query_proteins": selected_proteins,
                "score_threshold": score_threshold,
                "network_type": network_type,
            }
        except requests.exceptions.RequestException as exc:
            st.error(
                "無法連線 STRING API。請確認電腦目前可連上網路，"
                f"稍後再試。詳細訊息：{exc}"
            )
            st.stop()
        except Exception as exc:
            st.error(f"取得 STRING 資料時發生問題：{exc}")
            st.stop()

result = st.session_state["rbp_network_result"]

# Re-query when controls change, rather than accidentally showing stale results.
settings_match = (
    result["query_proteins"] == selected_proteins
    and result["score_threshold"] == score_threshold
    and result["network_type"] == network_type
)
if not settings_match:
    st.warning("查詢條件已變更。請按左側「查詢 STRING 網絡」更新結果。")
    st.stop()

mapping = result["mapping"]
network = result["network"]
missing = result["missing"]

if missing:
    st.warning("以下名稱未在 STRING 對應結果中找到：" + ", ".join(missing))

st.divider()
st.header("2. STRING protein ID mapping")

mapping_display = mapping.copy()
mapping_display["Biological context"] = mapping_display["preferredName"].map(
    PROTEIN_NOTES
).fillna("候選蛋白質")
st.dataframe(
    mapping_display[["queryItem", "preferredName", "stringId", "Biological context"]],
    use_container_width=True,
    hide_index=True,
)
st.caption(
    "STRING 會將輸入名稱對應至其內部 STRING protein identifier；"
    "請以回傳的 preferred name 與 STRING ID 作為後續網絡查詢依據。"
)

st.divider()
st.header("3. Protein network")

if network.empty:
    st.warning(
        "在目前候選蛋白、網絡類型及分數門檻下，STRING 沒有回傳連線。"
        "這不代表蛋白質之間絕對沒有關係；可嘗試降低分數門檻或改查 Functional association。"
    )
else:
    # STRING functional/physical TSV returns channel-specific evidence scores.
    expected = [
        "preferredName_A", "preferredName_B", "score",
        "nscore", "fscore", "pscore", "ascore",
        "escore", "dscore", "tscore"
    ]
    for col in expected:
        if col not in network.columns:
            network[col] = pd.NA

    # STRING network TSV commonly returns scores on a 0–1 scale.
    # Keep compatibility with responses that may return a 0–1000 scale.
    raw_score = pd.to_numeric(network["score"], errors="coerce")
    score_max = raw_score.max(skipna=True)
    if pd.notna(score_max) and score_max > 1:
        network["combined_score_0_1"] = raw_score / 1000.0
    else:
        network["combined_score_0_1"] = raw_score

    c1, c2, c3 = st.columns(3)
    c1.metric("Protein nodes", len(set(
        network["preferredName_A"].dropna().tolist()
        + network["preferredName_B"].dropna().tolist()
    )))
    c2.metric("Returned connections", len(network))
    c3.metric("Highest combined score", f"{network['combined_score_0_1'].max():.3f}")

    if nx is not None:
        graph = nx.Graph()
        for _, row in network.iterrows():
            a = str(row["preferredName_A"])
            b = str(row["preferredName_B"])
            score = float(row["combined_score_0_1"])
            graph.add_edge(a, b, weight=score)

        if graph.number_of_edges() > 0:
            fig, ax = plt.subplots(figsize=(11, 8))
            pos = nx.spring_layout(graph, seed=27, k=0.9)
            edge_weights = [
                0.8 + 3.0 * graph[u][v]["weight"]
                for u, v in graph.edges()
            ]
            nx.draw_networkx_edges(
                graph, pos, ax=ax, width=edge_weights, alpha=0.65
            )
            nx.draw_networkx_nodes(
                graph, pos, ax=ax, node_size=1800
            )
            nx.draw_networkx_labels(
                graph, pos, ax=ax, font_size=10, font_weight="bold"
            )
            ax.set_title(
                f"TARBP2 candidate network ({network_type_label})",
                fontsize=14,
            )
            ax.axis("off")
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)
    else:
        st.info(
            "目前環境未安裝 networkx，因此先顯示互作資料表。"
            "若要繪製網絡圖，可在終端機執行：python -m pip install networkx"
        )

    st.subheader("Interaction evidence table")
    display = network[
        [
            "preferredName_A", "preferredName_B", "score",
            "escore", "dscore", "tscore", "ascore",
            "nscore", "fscore", "pscore"
        ]
    ].copy()
    display = display.rename(
        columns={
            "preferredName_A": "Protein A",
            "preferredName_B": "Protein B",
            "score": "Combined score (0–1)",
            "escore": "Experimental",
            "dscore": "Database",
            "tscore": "Text mining",
            "ascore": "Co-expression",
            "nscore": "Neighborhood",
            "fscore": "Gene fusion",
            "pscore": "Phylogenetic profile",
        }
    )
    st.dataframe(
        display.sort_values(
            "Combined score (0–1)", ascending=False
        ),
        use_container_width=True,
        hide_index=True,
    )

    csv_bytes = network.to_csv(index=False).encode("utf-8-sig")
    st.download_button(
        "下載 STRING network CSV",
        data=csv_bytes,
        file_name="TARBP2_STRING_network.csv",
        mime="text/csv",
    )

st.divider()
st.header("4. How to interpret the evidence")

st.markdown(
    """
    | Evidence channel | 意義 |
    |---|---|
    | Combined score | STRING 整合各證據通道後的整體信心分數（0–1000） |
    | Experimental | 實驗資料支持的關聯證據 |
    | Database | 已整理資料庫中的路徑或互作證據 |
    | Text mining | 文獻文字探勘所支持的關聯 |
    | Co-expression | 基因共同表現的支持程度 |
    | Neighborhood / Gene fusion / Phylogenetic profile | 基因組背景或演化資訊所提供的支持 |

    **重要限制：**
    - STRING 的 combined score 是信心分數，不是結合親和力、相關係數或因果機率。
    - Functional association 不一定代表兩個蛋白質直接接觸。
    - Text mining 代表文獻文字中存在關聯線索，不等於已由實驗證實。
    - 本頁候選蛋白是依 TARBP2 / RNAi 生物學背景選定；網絡結果不能單獨證明 TARBP2 在特定癌症中的作用機制。
    """
)

st.markdown(
    """
    **Data source:** STRING database, Homo sapiens (NCBI taxonomy ID: 9606).
    API: identifier mapping and network retrieval.
    """
)
st.markdown(
    "官方資料："
    "[STRING API documentation](https://string-db.org/help/api/) · "
    "[STRING database access and licensing](https://string-db.org/cgi/access.pl)"
)
st.caption(
    "STRING data are made available under CC BY 4.0. "
    "When reusing outputs, cite STRING and indicate any changes or additions. "
    "This page retrieves data on demand and does not imply independent experimental validation."
)
