import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root / "build" / "python"))
import aocs_orbit

df = pd.read_csv("output.csv")
position = df[["r_x_m", "r_y_m", "r_z_m"]].to_numpy()
velocity = df[["v_x_m_s", "v_y_m_s", "v_z_m_s"]].to_numpy()
elements = aocs_orbit.convert_ijk_to_keplerian_series(position, velocity)
df["sma_m"] = elements[:, 0]
df["eccentricity"] = elements[:, 1]
df["inclination_deg"] = np.degrees(elements[:, 2])
df["raan_deg"] = np.degrees(elements[:, 3])
df["argp_deg"] = np.degrees(elements[:, 4])
df["true_anomaly_deg"] = np.degrees(elements[:, 5])

elements = [
    ("sma_m", "Semi-major axis", "m"),
    ("eccentricity", "Eccentricity", ""),
    ("inclination_deg", "Inclination", "deg"),
    ("raan_deg", "RAAN", "deg"),
    ("argp_deg", "Argument of periapsis", "deg"),
    ("true_anomaly_deg", "True anomaly", "deg"),
]

background = "#080b14"
panel = "#0d1422"
grid = "#1b2b40"
text = "#c5d5e5"
pearl = "#aeb1b0"
orange = "#fa8d11"
amber = "#ffd166"

fig = make_subplots(
    rows=len(elements) + 1,
    cols=1,
    specs=[[{"type": "scene"}]] + [[{"type": "xy"}] for _ in elements],
    row_heights=[0.48] + [0.52 / len(elements)] * len(elements),
    vertical_spacing=0.035,
    subplot_titles=["ECI orbit", *[title for _, title, _ in elements]],
)

earth_radius_km = 6371.0
latitude = np.linspace(-np.pi / 2.0, np.pi / 2.0, 48)
longitude = np.linspace(0.0, 2.0 * np.pi, 96)
latitude_grid, longitude_grid = np.meshgrid(latitude, longitude)
earth_x = earth_radius_km * np.cos(latitude_grid) * np.cos(longitude_grid)
earth_y = earth_radius_km * np.cos(latitude_grid) * np.sin(longitude_grid)
earth_z = earth_radius_km * np.sin(latitude_grid)

fig.add_trace(
    go.Surface(
        x=earth_x,
        y=earth_y,
        z=earth_z,
        surfacecolor=np.broadcast_to(latitude_grid, earth_x.shape),
        colorscale=[[0.0, "#10243a"], [0.5, "#22526E"], [1.0, "#247c83"]],
        cmin=-np.pi / 2.0,
        cmax=np.pi / 2.0,
        showscale=False,
        opacity=0.98,
        lighting={
            "ambient": 0.2,
            "diffuse": 0.9,
            "specular": 0.25,
            "roughness": 0.75,
            "fresnel": 0.0125,
        },
        lightposition={"x": 10000, "y": 12000, "z": 18000},
        hovertemplate="Earth<extra></extra>",
        name="Earth",
    ),
    row=1,
    col=1,
)

position_km = position / 1000.0
fig.add_trace(
    go.Scatter3d(
        x=position_km[:, 0],
        y=position_km[:, 1],
        z=position_km[:, 2],
        mode="lines",
        line={"color": pearl, "width": 4},
        name="Orbit",
        hovertemplate="x=%{x:.1f} km<br>y=%{y:.1f} km<br>z=%{z:.1f} km<extra>Orbit</extra>",
    ),
    row=1,
    col=1,
)
fig.add_trace(
    go.Scatter3d(
        x=[position_km[0, 0]],
        y=[position_km[0, 1]],
        z=[position_km[0, 2]],
        mode="markers",
        marker={"color": amber, "size": 5, "symbol": "circle"},
        name="Start",
        hovertemplate="Start<extra></extra>",
    ),
    row=1,
    col=1,
)
fig.add_trace(
    go.Scatter3d(
        x=[position_km[-1, 0]],
        y=[position_km[-1, 1]],
        z=[position_km[-1, 2]],
        mode="markers",
        marker={"color": orange, "size": 6, "symbol": "diamond"},
        name="End",
        hovertemplate="End<extra></extra>",
    ),
    row=1,
    col=1,
)

for row, (column, title, unit) in enumerate(elements, start=2):
    fig.add_trace(
        go.Scatter(
            x=df["t_sec"],
            y=df[column],
            mode="lines",
            name=title,
            line={"color": pearl, "width": 2},
            hovertemplate="%{x:.1f} s<br>%{y:.5g} " + unit + "<extra></extra>",
        ),
        row=row,
        col=1,
    )
    fig.update_yaxes(
        title_text=unit,
        row=row,
        col=1,
        gridcolor=grid,
        zerolinecolor=grid,
        linecolor=grid,
        tickfont={"color": text},
        title_font={"color": text},
    )

fig.update_xaxes(
    title_text="Time (s)",
    row=len(elements) + 1,
    col=1,
    gridcolor=grid,
    zerolinecolor=grid,
    linecolor=grid,
    tickfont={"color": text},
    title_font={"color": text},
)
fig.update_layout(
    title={
        "text": "ORBITAL TELEMETRY  //  OSCULATING ELEMENTS",
        "x": 0.04,
        "font": {"size": 22, "color": pearl, "family": "monospace"},
    },
    template="plotly_dark",
    paper_bgcolor=background,
    plot_bgcolor=panel,
    font={"color": text, "family": "monospace"},
    height=1850,
    margin={"l": 80, "r": 45, "t": 95, "b": 55},
    showlegend=False,
    hoverlabel={"bgcolor": panel, "bordercolor": pearl, "font": {"color": text}},
    scene={
        "bgcolor": background,
        "xaxis": {
            "title": "X (km)",
            "backgroundcolor": panel,
            "gridcolor": grid,
            "zerolinecolor": grid,
            "color": text,
            "showbackground": True,
        },
        "yaxis": {
            "title": "Y (km)",
            "backgroundcolor": panel,
            "gridcolor": grid,
            "zerolinecolor": grid,
            "color": text,
            "showbackground": True,
        },
        "zaxis": {
            "title": "Z (km)",
            "backgroundcolor": panel,
            "gridcolor": grid,
            "zerolinecolor": grid,
            "color": text,
            "showbackground": True,
        },
        "aspectmode": "data",
        "camera": {"eye": {"x": 1.55, "y": 1.55, "z": 1.2}},
    },
)
fig.update_annotations(font={"color": orange, "size": 13, "family": "monospace"})
fig.write_html("output.html")
