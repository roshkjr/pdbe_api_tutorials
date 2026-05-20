import marimo

__generated_with = "0.23.6"
app = marimo.App(width="medium")


@app.cell
def _():
    import altair as alt
    import marimo as mo
    import math
    import pandas as pd
    import requests
    import textwrap

    return alt, math, mo, pd, requests, textwrap


@app.cell
def _(mo):
    mo.md(r"""
    [![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/roshkjr/pdbe_api_tutorials/blob/main/notebooks/egfr_ligand_batch_analysis.py)
    """)
    return


@app.cell
def _(mo, textwrap):
    mo.md(
        textwrap.dedent(
            """
            # EGFR ligand analysis with PDBe APIs

            This notebook follows analysis of ligands in therapeutic target **EGFR**
            (`UniProt: P00533`):

            1. Fetch the ligands observed for EGFR with a **GET** request to the `/uniprot/ligands/` endpoint.
            2. Use a **POST** request to fetch compound summaries in batch from `/pdb/compound/summary` endpoint.
            3. Keep only ligands with a **DrugBank** cross-reference and number of heavy atoms greater than 6
            4. Summarize how often those ligands appear across EGFR structures.
            5. Group the retained ligands by Murcko scaffolds.

            API docs: <https://www.ebi.ac.uk/pdbe/api/v2/doc/>
            """
        )
    )
    return


@app.cell
def _():
    BASE_URL = "https://www.ebi.ac.uk/pdbe/api/v2"
    DOCS_URL = "https://www.ebi.ac.uk/pdbe/api/v2/doc/"
    TARGET_NAME = "EGFR"
    TARGET_UNIPROT_ID = "P00533"
    return BASE_URL, DOCS_URL, TARGET_NAME, TARGET_UNIPROT_ID


@app.cell
def _(BASE_URL, requests):
    session = requests.Session()
    session.headers.update(
        {
            "Accept": "application/json",
            "User-Agent": "pdbe-api-tutorials-marimo/0.1.0",
        }
    )

    def request_json(method, path, *, json_payload=None, timeout=60):
        url = f"{BASE_URL}{path}"
        response = session.request(
            method=method,
            url=url,
            json=json_payload,
            timeout=timeout,
        )
        response.raise_for_status()
        return response, response.json()

    def get_json(path, *, timeout=60):
        return request_json("GET", path, timeout=timeout)

    def post_json(path, payload, *, timeout=60):
        return request_json("POST", path, json_payload=payload, timeout=timeout)


    def summarize_response(response, *, label):
        return {
            "request": label,
            "status_code": response.status_code,
            "content_type": response.headers.get("content-type", ""),
            "elapsed_ms": round(response.elapsed.total_seconds() * 1000, 1),
            "url": response.url,
        }

    return get_json, post_json, session, summarize_response


@app.cell
def _(TARGET_NAME, TARGET_UNIPROT_ID, mo):
    mo.md(f"""
    ## 1. Fetch EGFR ligands with GET

    The first live request starts from the biological target:

    - Endpoint: `/uniprot/ligands/{TARGET_UNIPROT_ID}`
    - Method: `GET`
    - Meaning: ligands observed for **{TARGET_NAME}** in PDBe, together with the PDB entries
      where they appear.
    """)
    return


@app.cell
def _(TARGET_UNIPROT_ID, get_json, summarize_response):
    uniprot_ligands_response, uniprot_ligands_json = get_json(
        f"/uniprot/ligands/{TARGET_UNIPROT_ID}"
    )
    uniprot_ligands_summary = summarize_response(
        uniprot_ligands_response,
        label="GET /uniprot/ligands/{uniprot_accession}",
    )
    return uniprot_ligands_json, uniprot_ligands_summary


@app.cell
def _(pd, uniprot_ligands_summary):
    pd.DataFrame([uniprot_ligands_summary])
    return


@app.cell
def _(TARGET_UNIPROT_ID, uniprot_ligands_json):
    uniprot_ligands_json[TARGET_UNIPROT_ID][0]
    return


@app.cell
def _(TARGET_UNIPROT_ID, pd, uniprot_ligands_json):
    ligand_rows = []
    for ligand_item in uniprot_ligands_json.get(TARGET_UNIPROT_ID, []):
        _hetcode, _payload = next(iter(ligand_item.items()))
        pdb_ids = sorted(set(_payload.get("pdbs", [])))
        acts_as = _payload.get("acts_as", [])
        ligand_rows.append(
            {
                "ligand_id": _hetcode,
                "ligand_name": _payload.get("name"),
                "pdb_ids": pdb_ids,
                "pdb_entry_count": len(pdb_ids),
                "acts_as": ", ".join(acts_as) if acts_as else "not annotated",
                "scaffold_id": _payload.get("scaffold_id"),
                "fragment_count": len(_payload.get("fragments", [])),
                "directly_interacts": _payload.get("directly_interacts"),
            }
        )

    ligands_df = pd.DataFrame(ligand_rows).sort_values(
        ["pdb_entry_count"],
        ascending=[True],
    )
    unique_ligand_ids = ligands_df["ligand_id"].tolist()
    unique_pdb_ids = sorted(
        {pdb_id for pdb_ids in ligands_df["pdb_ids"] for pdb_id in pdb_ids}
    )
    return ligands_df, unique_ligand_ids, unique_pdb_ids


@app.cell
def _(ligands_df):
    ligands_df.head(15)
    return


@app.cell
def _(TARGET_NAME, ligands_df, mo, unique_pdb_ids):
    mo.md(f"""
    The live GET response currently contains:

    - **{len(ligands_df):,}** unique ligand identifiers for {TARGET_NAME}
    - **{len(unique_pdb_ids):,}** unique PDB entries represented in those ligand records
    """)
    return


@app.cell
def _(mo, textwrap):
    mo.md(
        textwrap.dedent(
            """
            ## 2. Why a large GET is brittle

            The compound summary API has a convenient batch **POST** route, but a comma-separated
            ligand list also works in the path for **GET**. Let's see when it won't work.
            """
        )
    )
    return


@app.cell
def _(BASE_URL, math, pd, session, unique_ligand_ids):
    demo_hetcodes = unique_ligand_ids[:20] if len(unique_ligand_ids) >= 20 else unique_ligand_ids
    if not demo_hetcodes:
        demo_hetcodes = ["ATP"]

    demo_joined = ",".join(demo_hetcodes)
    target_url_length = 20000
    multiplier = max(1, math.ceil(target_url_length / max(len(demo_joined), 1)))
    oversized_hetcodes = demo_hetcodes * multiplier
    oversized_get_url = f"{BASE_URL}/pdb/compound/summary/{','.join(oversized_hetcodes)}"

    oversized_get_response = session.get(oversized_get_url, timeout=60)
    oversized_get_result = {
        "status_code": oversized_get_response.status_code,
        "url_length": len(oversized_get_url),
        "ligand_tokens_in_url": len(oversized_hetcodes),
        "body_snippet": oversized_get_response.text[:250],
    }

    pd.DataFrame([oversized_get_result]).drop(columns=["body_snippet"])
    return (oversized_get_result,)


@app.cell
def _(mo, oversized_get_result):
    mo.md(f"""
    Oversized GET body preview:

    ```text
    {oversized_get_result["body_snippet"]}
    ```

    In practice, this is the point of the batch **POST** route: once the identifier set grows,
    the URL becomes fragile and can fail with **`414 Request-URI Too Large`**.
    """)
    return


@app.cell
def _(mo, textwrap):
    mo.md(
        textwrap.dedent(
            """
            ## 3. Batch fetch compound summaries with POST

            The PDBe batch summary route is:

            - Endpoint: `/pdb/compound/summary`
            - Method: `POST`
            - Request body: a JSON **string** containing comma-separated hetcodes

            This keeps the URL short while moving the identifier list into the request body.
            """
        )
    )
    return


@app.cell
def _(mo, unique_ligand_ids):
    preview_payload = ",".join(unique_ligand_ids[:12])
    mo.md(
        f"""
        Example POST payload shape:

        ```json
        "{preview_payload}"
        ```
        """
    )
    return


@app.cell
def _(post_json, summarize_response, unique_ligand_ids):
    compound_summary_payload = ",".join(unique_ligand_ids)
    compound_summary_response, compound_summary_json = post_json(
        "/pdb/compound/summary",
        compound_summary_payload,
    )
    compound_summary_summary = summarize_response(
        compound_summary_response,
        label="POST /pdb/compound/summary",
    )
    return compound_summary_json, compound_summary_summary


@app.cell
def _(compound_summary_summary, pd):
    pd.DataFrame([compound_summary_summary])
    return


@app.cell
def _(compound_summary_json):
    first_summary_keys = list(compound_summary_json.keys())[:1]
    raw_compound_preview = {
        key: compound_summary_json[key] for key in first_summary_keys
    }
    raw_compound_preview
    return


@app.cell
def _(compound_summary_json, ligands_df, pd):
    summary_rows = []
    for _, ligand_row in ligands_df.iterrows():
        _summary_hetcode = ligand_row["ligand_id"]
        summary_payload = compound_summary_json.get(_summary_hetcode, [])
        summary = summary_payload[0] if summary_payload else {}
        cross_links = summary.get("cross_links") or []
        smiles_payload = summary.get("smiles") or []
        primary_smiles = next(
            (entry.get("name") for entry in smiles_payload if entry.get("name")),
            None,
        )
        drugbank_ids = sorted(
            {
                link.get("resource_id")
                for link in cross_links
                if link.get("resource") == "DrugBank" and link.get("resource_id")
            }
        )
        drugbank_synonyms = [
            synonym.get("value")
            for synonym in (summary.get("synonyms") or [])
            if synonym.get("origin") == "DrugBank" and synonym.get("value")
        ]
        phys_chem = summary.get("phys_chem_properties") or {}

        summary_rows.append(
            {
                "ligand_id": _summary_hetcode,
                "summary_found": bool(summary),
                "compound_name": summary.get("name") or ligand_row["ligand_name"],
                "weight": summary.get("weight"),
                "exact_mw": phys_chem.get("exactmw"),
                "tpsa": phys_chem.get("tpsa"),
                "clogp": phys_chem.get("crippen_clog_p"),
                "smiles": primary_smiles,
                "heavy_atoms": phys_chem.get("num_heavy_atoms"),
                "hba": phys_chem.get("num_hba"),
                "hbd": phys_chem.get("num_hbd"),
                "drugbank_ids": ", ".join(drugbank_ids),
                "drugbank_id_count": len(drugbank_ids),
                "drugbank_synonyms": ", ".join(drugbank_synonyms),
                "cross_link_count": len(cross_links),
            }
        )

    compound_summary_df = ligands_df.merge(
        pd.DataFrame(summary_rows),
        on="ligand_id",
        how="left",
    )
    return (compound_summary_df,)


@app.cell
def _(compound_summary_df, pd):
    counts_df = pd.DataFrame(
        [
            {
                "metric": "Unique ligands observed for EGFR",
                "value": int(compound_summary_df["ligand_id"].nunique()),
            },
            {
                "metric": "Ligands with compound summary records",
                "value": int(compound_summary_df["summary_found"].fillna(False).sum()),
            },
            {
                "metric": "Ligands retained after DrugBank filter",
                "value": int((compound_summary_df["drugbank_id_count"].fillna(0) > 0).sum()),
            },
        ]
    )
    counts_df
    return


@app.cell
def _(compound_summary_df):
    compound_summary_df
    return


@app.cell
def _(compound_summary_df):
    filtered_ligands_df = compound_summary_df.loc[
        (compound_summary_df["drugbank_id_count"].fillna(0) > 0)
        & (compound_summary_df["heavy_atoms"].fillna(0) > 6)
    ].copy()
    filtered_ligands_df["example_pdb_entries"] = filtered_ligands_df["pdb_ids"].apply(
        lambda pdb_ids: ", ".join(pdb_ids[:6]) if pdb_ids else ""
    )
    filtered_ligands_df = filtered_ligands_df.sort_values(
        ["pdb_entry_count", "weight", "ligand_id"],
        ascending=[False, False, True],
    )

    display_columns = [
        "ligand_id",
        "compound_name",
        "drugbank_ids",
        "pdb_entry_count",
        "weight",
        "acts_as",
        "example_pdb_entries",
    ]
    filtered_ligands_display_df = filtered_ligands_df[display_columns].reset_index(drop=True)
    return filtered_ligands_df, filtered_ligands_display_df


@app.cell
def _(filtered_ligands_df, mo):
    if filtered_ligands_df.empty:
        filtered_ligands_message = mo.md(
            """
            No EGFR ligands with a DrugBank cross-reference were returned in this live run.

            The earlier cells still show the fetched EGFR ligands, the POST summary enrichment,
            and the pre/post-filter counts. The archive and cross-references are live and may
            change over time.
            """
        )
    else:
        filtered_ligands_message = mo.md(
            """
            ## 4. DrugBank-filtered ligands

            The table below keeps only ligands whose PDBe compound summary includes at least one
            **DrugBank** cross-reference and at least 7 heavy atoms.
            """
        )
    filtered_ligands_message
    return


@app.cell
def _(filtered_ligands_display_df):
    filtered_ligands_display_df
    return


@app.cell
def _(mo, textwrap):
    mo.md(
        textwrap.dedent(
            """
            ## 5. Analysis 1: recurrence across EGFR structures

            A simple first analysis is to count how many EGFR-associated PDB entries contain each
            retained ligand.
            """
        )
    )
    return


@app.cell
def _(alt, filtered_ligands_df, mo):
    if filtered_ligands_df.empty:
        ligand_frequency_chart = None
        ligand_frequency_view = mo.md(
            "No bar chart to render because no DrugBank-mapped ligands were retained."
        )
    else:
        chart_df = filtered_ligands_df.nlargest(15, "pdb_entry_count")
        ligand_frequency_chart = (
            alt.Chart(chart_df)
            .mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4, color="#1f77b4")
            .encode(
                x=alt.X(
                    "pdb_entry_count:Q",
                    title="Number of EGFR-associated PDB entries",
                ),
                y=alt.Y(
                    "ligand_id:N",
                    sort="-x",
                    title="Ligand",
                ),
                tooltip=[
                    alt.Tooltip("compound_name:N", title="Ligand"),
                    alt.Tooltip("drugbank_ids:N", title="DrugBank"),
                    alt.Tooltip("pdb_entry_count:Q", title="PDB entries"),
                    alt.Tooltip("weight:Q", title="Formula weight"),
                ],
            )
            .properties(
                title="DrugBank-mapped ligands observed across EGFR structures",
                height=420,
            )
        )
        ligand_frequency_view = ligand_frequency_chart
    ligand_frequency_view
    return


@app.cell
def _(mo, textwrap):
    mo.md(
        textwrap.dedent(
            """
            ## 6. Analysis 2: scaffold grouping from the API

            The UniProt ligand response includes a precomputed `scaffold_id` field. This corresponds
            to Murcko scaffold calcualted using RDKit. This section groups the retained ligands by their scaffolds.
            """
        )
    )
    return


@app.cell
def _(filtered_ligands_df):
    scaffolded_ligands_df = filtered_ligands_df.copy()
    scaffolded_ligands_df["murcko_scaffold_group"] = scaffolded_ligands_df["scaffold_id"].fillna(
        "No scaffold assigned"
    )
    scaffolded_ligands_df["murcko_scaffold_member_count"] = scaffolded_ligands_df.groupby(
        "murcko_scaffold_group"
    )["ligand_id"].transform("nunique")
    scaffolded_ligands_df["murcko_scaffold_total_pdb_entries"] = scaffolded_ligands_df.groupby(
        "murcko_scaffold_group"
    )["pdb_entry_count"].transform("sum")

    scaffold_group_df = (
        scaffolded_ligands_df.groupby("murcko_scaffold_group", dropna=False)
        .agg(
            ligand_count=("ligand_id", "nunique"),
            total_pdb_entry_count=("pdb_entry_count", "sum"),
            max_single_ligand_pdb_count=("pdb_entry_count", "max"),
            ligand_ids=("ligand_id", lambda values: ", ".join(sorted(values))),
            ligand_names=("compound_name", lambda values: ", ".join(sorted(set(values)))),
        )
        .reset_index()
        .sort_values(
            ["ligand_count", "total_pdb_entry_count", "murcko_scaffold_group"],
            ascending=[False, False, True],
        )
    )
    return scaffold_group_df, scaffolded_ligands_df


@app.cell
def _(alt, mo, scaffold_group_df):
    if scaffold_group_df.empty:
        scaffold_chart = None
        scaffold_chart_view = mo.md(
            "No scaffold grouping chart to render because no DrugBank-mapped ligands were retained."
        )
    else:
        scaffold_chart = (
            alt.Chart(scaffold_group_df.head(15))
            .mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4, color="#0f766e")
            .encode(
                x=alt.X("ligand_count:Q", title="Number of ligands in scaffold group"),
                y=alt.Y("murcko_scaffold_group:N", sort="-x", title="Murcko scaffold"),
                tooltip=[
                    alt.Tooltip("murcko_scaffold_group:N", title="Scaffold"),
                    alt.Tooltip("ligand_count:Q", title="Ligands"),
                    alt.Tooltip("total_pdb_entry_count:Q", title="Total PDB entries"),
                    alt.Tooltip("ligand_ids:N", title="Ligand IDs"),
                ],
            )
            .properties(
                title="Murcko scaffold groups among DrugBank-mapped EGFR ligands",
                height=420,
            )
        )
        scaffold_chart_view = scaffold_chart
    scaffold_chart_view
    return


@app.cell
def _(scaffold_group_df):
    scaffold_group_df.reset_index(drop=True)
    return


@app.cell
def _(scaffolded_ligands_df):
    scaffold_member_columns = [
        "ligand_id",
        "compound_name",
        "drugbank_ids",
        "pdb_entry_count",
        "murcko_scaffold_group",
        "murcko_scaffold_member_count",
        "murcko_scaffold_total_pdb_entries",
    ]
    scaffolded_ligands_df[scaffold_member_columns].sort_values(
        ["murcko_scaffold_member_count", "pdb_entry_count", "ligand_id"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    return


@app.cell
def _(DOCS_URL, TARGET_NAME, TARGET_UNIPROT_ID, filtered_ligands_df, mo):
    retained_count = len(filtered_ligands_df)
    mo.md(
        f"""
        ## Takeaways

        - The biological entry point is a target-centric **GET** request:
          `/uniprot/ligands/{TARGET_UNIPROT_ID}` for **{TARGET_NAME}**.
        - Trying to batch too many ligand identifiers into a **GET** path can fail with
          **`414 Request-URI Too Large`**.
        - The PDBe batch route `POST /pdb/compound/summary` is the robust pattern for enriching
          many ligands at once.
        - The UniProt ligand payload already provides a scaffold field, which can be used directly
          to group retained ligands by shared core chemistry.
        - In this live run, **{retained_count}** EGFR-associated ligands were retained after
          requiring a **DrugBank** cross-reference.

        Live docs: {DOCS_URL}
        """
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ## Try yourself

    1. Find out which residues are mostly interacting with EGFR ligands across structures ?
    2. What type of interactions are mostly found ?
    3. Are there any mutated residues in the ligand binding sites?
    """)
    return


@app.cell
def _(mo):
    mo.vstack([mo.md("## Hints"),
    mo.accordion({
        "Get all the ligands bound to an entry": "`/pdb/bound_molecules/<pdb_id>`",
        "For ligand interactions in an entry": "`/pdb/bound_ligand_interactions/<pdb_id>/<chain_id>/<seq_id>`",
        "For mutations see the endpoint": "`/uniprot/unipdb/<uniprot_id>`",
        "For mapping between uniprot residues and PDB residues": "`/mappings/uniprot/<pdb_id>`",
        }

    )
    ])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
