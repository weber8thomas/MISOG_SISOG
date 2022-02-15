from dash import Dash, html, dcc, Input, Output, dash_table
import pandas as pd
import plotly.express as px
import dash_bio

external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

app = Dash(__name__, external_stylesheets=external_stylesheets)


fig2 = dash_bio.Igv()

# df = pd.read_parquet("/gstock/GeneIso/V2/Introns.parquet")
# df["Length"] = df["Introns"].apply(lambda r: int(r.split("-")[1]) - int(r.split("-")[0]))
# df = df.loc[df["Introns_nb"] < 6]
# df = pd.concat([df.loc[df["Miso_siso"] == "Siso"].head(500), df.loc[df["Miso_siso"] == "Miso"].head(500)])
# print(df)


app.layout = html.Div(
    [
        html.H1("TEST"),
        dash_bio.Igv(
            id="reference-igv",
            reference={
                "id": "ASM985889v3",
                "name": "Sars-CoV-2 (ASM985889v3)",
                "fastaURL": "https://s3.amazonaws.com/igv.org.genomes/covid_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna",
                "indexURL": "https://s3.amazonaws.com/igv.org.genomes/covid_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.fai",
                "order": 1000000,
                "tracks": [
                    {
                        "name": "Annotations",
                        "url": "https://s3.amazonaws.com/igv.org.genomes/covid_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz",
                        "displayMode": "EXPANDED",
                        "nameField": "gene",
                        "height": 150,
                        "color": "rgb(176,141,87)",
                    }
                ],
            },
        ),
    ]
)
#         html.H1("TEST"),
#         html.Div(
#             [
#                 dcc.Dropdown(df.columns.tolist(), id="xaxis", value="Introns_nb"),
#                 dcc.Dropdown(df.columns.tolist(), id="yaxis", value="Length"),
#                 dcc.Dropdown(df.columns.tolist(), id="color", value="Miso_siso"),
#             ]
#         ),
#         dcc.Graph(id="graph"),
#         html.Div(id="table"),
#     ]
# )


# @app.callback(
#     Output("graph", "figure"),
#     Input("xaxis", "value"),
#     Input("yaxis", "value"),
#     Input("color", "value"),
# )
# def figure(xaxis, yaxis, color):
#     return px.violin(df, x=xaxis, y=yaxis, color=color)


# @app.callback(
#     Output("table", "children"),
#     Input("color", "value"),
# )
# def table(x):
#     return dash_table.DataTable(df.head(20).to_dict("records"))


if __name__ == "__main__":
    app.run_server(debug=True)
