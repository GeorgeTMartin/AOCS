import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root / "build" / "python"))
import aocs_orbit


BACKGROUND = "#080b14"
PANEL = "#0d1422"
GRID = "#1b2b40"
TEXT = "#c5d5e5"
PEARL = "#aeb1b0"
ORANGE = "#fa8d11"
AMBER = "#ffd166"
CYAN = "#48d6d2"
GREEN = "#9ce86b"
BLUE = "#58a6ff"
PLOTLY_CONFIG = {"responsive": True, "displaylogo": False}


def build_state_figure(df):
    position = df[["r_x_m", "r_y_m", "r_z_m"]].to_numpy()
    velocity = df[["v_x_m_s", "v_y_m_s", "v_z_m_s"]].to_numpy()
    elements = aocs_orbit.convert_ijk_to_keplerian_series(position, velocity)
    df["sma_m"] = elements[:, 0]
    df["eccentricity"] = elements[:, 1]
    df["inclination_deg"] = np.degrees(elements[:, 2])
    df["raan_deg"] = np.degrees(elements[:, 3])
    df["argp_deg"] = np.degrees(elements[:, 4])
    df["true_anomaly_deg"] = np.degrees(elements[:, 5])

    specs = [
        [{"type": "scene", "rowspan": 3}, {"type": "xy"}, {"type": "xy"}],
        [None, {"type": "xy"}, {"type": "xy"}],
        [None, {"type": "xy"}, {"type": "xy"}],
        [{"type": "xy"}, {"type": "xy"}, {"type": "xy"}],
        [{"type": "xy"}, {"type": "xy"}, {"type": "xy"}],
    ]
    titles = [
        "ECI Orbit",
        "X Position",
        "X Velocity",
        "",
        "Y Position",
        "Y Velocity",
        "",
        "Z Position",
        "Z Velocity",
        "Semi-major Axis",
        "Eccentricity",
        "Inclination",
        "RAAN",
        "Argument of Periapsis",
        "True Anomaly",
    ]
    fig = make_subplots(
        rows=5,
        cols=3,
        specs=specs,
        row_heights=[0.18, 0.18, 0.18, 0.23, 0.23],
        column_widths=[0.40, 0.30, 0.30],
        horizontal_spacing=0.055,
        vertical_spacing=0.055,
        subplot_titles=titles,
    )

    latitude = np.linspace(-np.pi / 2.0, np.pi / 2.0, 40)
    longitude = np.linspace(0.0, 2.0 * np.pi, 80)
    latitude_grid, longitude_grid = np.meshgrid(latitude, longitude)
    earth_radius_km = 6371.0
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
            line={"color": PEARL, "width": 4},
            name="Orbit",
            hovertemplate=(
                "x=%{x:.1f} km<br>y=%{y:.1f} km<br>"
                "z=%{z:.1f} km<extra>Orbit</extra>"
            ),
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
            marker={"color": AMBER, "size": 5},
            name="Start",
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
            marker={"color": ORANGE, "size": 6, "symbol": "diamond"},
            name="End",
        ),
        row=1,
        col=1,
    )

    position_colors = [CYAN, GREEN, BLUE]
    velocity_colors = [ORANGE, AMBER, PEARL]
    position_names = ["X", "Y", "Z"]
    velocity_names = ["Vx", "Vy", "Vz"]
    for axis in range(3):
        fig.add_trace(
            go.Scatter(
                x=df["t_sec"],
                y=position[:, axis] / 1000.0,
                mode="lines",
                line={"color": position_colors[axis], "width": 1.5},
                name=f"{position_names[axis]} position",
                hovertemplate="%{x:.1f} s<br>%{y:.4f} km<extra></extra>",
            ),
            row=axis + 1,
            col=2,
        )
        fig.add_trace(
            go.Scatter(
                x=df["t_sec"],
                y=velocity[:, axis] / 1000.0,
                mode="lines",
                line={"color": velocity_colors[axis], "width": 1.5},
                name=f"{velocity_names[axis]} velocity",
                hovertemplate="%{x:.1f} s<br>%{y:.5f} km/s<extra></extra>",
            ),
            row=axis + 1,
            col=3,
        )
        fig.update_yaxes(
            title_text="km",
            row=axis + 1,
            col=2,
            gridcolor=GRID,
            zerolinecolor=GRID,
            linecolor=GRID,
            tickfont={"color": TEXT},
            title_font={"color": TEXT},
        )
        fig.update_yaxes(
            title_text="km/s",
            row=axis + 1,
            col=3,
            gridcolor=GRID,
            zerolinecolor=GRID,
            linecolor=GRID,
            tickfont={"color": TEXT},
            title_font={"color": TEXT},
        )

    element_panels = [
        ("sma_m", "Semi-major axis", "km", lambda values: values / 1000.0),
        ("eccentricity", "Eccentricity", "", lambda values: values),
        ("inclination_deg", "Inclination", "deg", lambda values: values),
        ("raan_deg", "RAAN", "deg", lambda values: values),
        ("argp_deg", "Argument of periapsis", "deg", lambda values: values),
        ("true_anomaly_deg", "True anomaly", "deg", lambda values: values),
    ]
    for index, (column, name, unit, conversion) in enumerate(element_panels):
        row = 4 + index // 3
        col = 1 + index % 3
        fig.add_trace(
            go.Scatter(
                x=df["t_sec"],
                y=conversion(df[column].to_numpy()),
                mode="lines",
                line={"color": PEARL, "width": 1.5},
                name=name,
                hovertemplate=f"%{{x:.1f}} s<br>%{{y:.5g}} {unit}<extra></extra>",
            ),
            row=row,
            col=col,
        )
        fig.update_yaxes(
            title_text=unit,
            row=row,
            col=col,
            gridcolor=GRID,
            zerolinecolor=GRID,
            linecolor=GRID,
            tickfont={"color": TEXT},
            title_font={"color": TEXT},
        )

    fig.update_xaxes(
        title_text="Time (s)",
        gridcolor=GRID,
        zerolinecolor=GRID,
        linecolor=GRID,
        tickfont={"color": TEXT},
        title_font={"color": TEXT},
    )
    fig.update_layout(
        title={
            "text": "ORBITAL TELEMETRY  //  STATE & OSCULATING ELEMENTS",
            "x": 0.035,
            "font": {"size": 22, "color": PEARL, "family": "monospace"},
        },
        template="plotly_dark",
        paper_bgcolor=BACKGROUND,
        plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "monospace"},
        height=1420,
        margin={"l": 65, "r": 35, "t": 95, "b": 55},
        showlegend=False,
        hoverlabel={"bgcolor": PANEL, "bordercolor": PEARL, "font": {"color": TEXT}},
        scene={
            "bgcolor": BACKGROUND,
            "xaxis": {"title": "X (km)", "backgroundcolor": PANEL, "gridcolor": GRID},
            "yaxis": {"title": "Y (km)", "backgroundcolor": PANEL, "gridcolor": GRID},
            "zaxis": {"title": "Z (km)", "backgroundcolor": PANEL, "gridcolor": GRID},
            "aspectmode": "data",
            "camera": {"eye": {"x": 1.55, "y": 1.55, "z": 1.2}},
        },
    )
    fig.update_annotations(font={"color": ORANGE, "size": 11, "family": "monospace"})
    return fig


def covariance_columns():
    return [f"cov_{row}_{col}" for row in range(6) for col in range(6)]


def ellipsoid_surface(covariance, scale, unit_x, unit_y, unit_z):
    symmetric_covariance = 0.5 * (covariance + covariance.T)
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric_covariance)
    radii = scale * np.sqrt(np.maximum(eigenvalues, 0.0))
    transform = eigenvectors @ np.diag(radii)
    unit_points = np.stack((unit_x.ravel(), unit_y.ravel(), unit_z.ravel()))
    points = transform @ unit_points
    return tuple(points[axis].reshape(unit_x.shape) for axis in range(3))


def build_covariance_figure(df):
    covariance = df[covariance_columns()].to_numpy().reshape(-1, 6, 6)
    if not np.isfinite(covariance).all():
        raise ValueError("Covariance data contains non-finite values")

    fig = make_subplots(
        rows=3,
        cols=1,
        specs=[[{"type": "xy"}], [{"type": "xy"}], [{"type": "scene"}]],
        row_heights=[0.19, 0.19, 0.62],
        vertical_spacing=0.08,
        subplot_titles=[
            "Position Variance",
            "Velocity Variance",
            "Animated 3σ Position Covariance Ellipsoid (local ECI offset)",
        ],
    )
    axis_colors = [CYAN, GREEN, BLUE]
    axis_names = ["X", "Y", "Z"]
    for axis in range(3):
        fig.add_trace(
            go.Scatter(
                x=df["t_sec"],
                y=covariance[:, axis, axis],
                mode="lines",
                line={"color": axis_colors[axis], "width": 1.5},
                name=f"P{axis_names[axis]}{axis_names[axis]}",
                hovertemplate="%{x:.1f} s<br>%{y:.6g} m²<extra></extra>",
            ),
            row=1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=df["t_sec"],
                y=covariance[:, axis + 3, axis + 3],
                mode="lines",
                line={"color": axis_colors[axis], "width": 1.5},
                name=f"V{axis_names[axis]}{axis_names[axis]}",
                hovertemplate="%{x:.1f} s<br>%{y:.6g} (m/s)²<extra></extra>",
            ),
            row=2,
            col=1,
        )
    for row, unit in ((1, "m²"), (2, "(m/s)²")):
        fig.update_xaxes(
            title_text="Time (s)",
            row=row,
            col=1,
            gridcolor=GRID,
            zerolinecolor=GRID,
            linecolor=GRID,
            tickfont={"color": TEXT},
            title_font={"color": TEXT},
        )
        fig.update_yaxes(
            title_text=unit,
            row=row,
            col=1,
            gridcolor=GRID,
            zerolinecolor=GRID,
            linecolor=GRID,
            tickfont={"color": TEXT},
            title_font={"color": TEXT},
        )

    latitude = np.linspace(0.0, np.pi, 18)
    longitude = np.linspace(0.0, 2.0 * np.pi, 36)
    latitude_grid, longitude_grid = np.meshgrid(latitude, longitude, indexing="ij")
    unit_x = np.sin(latitude_grid) * np.cos(longitude_grid)
    unit_y = np.sin(latitude_grid) * np.sin(longitude_grid)
    unit_z = np.cos(latitude_grid)
    frame_indices = np.unique(
        np.linspace(0, len(df) - 1, min(len(df), 120), dtype=int)
    )
    estimate_positions = df[
        ["estimate_r_x_m", "estimate_r_y_m", "estimate_r_z_m"]
    ].to_numpy()
    if not np.isfinite(estimate_positions).all():
        raise ValueError("Estimated position data contains non-finite values")

    frames = []
    initial_surface = None
    for frame_index in frame_indices:
        x, y, z = ellipsoid_surface(
            covariance[frame_index, :3, :3], 3.0, unit_x, unit_y, unit_z
        )
        if initial_surface is None:
            initial_surface = (x, y, z)
        time_s = df["t_sec"].iloc[frame_index]
        frames.append(
            go.Frame(
                name=f"{time_s:.2f}",
                data=[
                    go.Surface(
                        x=x,
                        y=y,
                        z=z,
                        colorscale=[[0.0, CYAN], [1.0, ORANGE]],
                        opacity=0.72,
                        showscale=False,
                        hovertemplate=(
                            "ΔX=%{x:.3f} m<br>ΔY=%{y:.3f} m"
                            "<br>ΔZ=%{z:.3f} m<extra>3σ ellipsoid</extra>"
                        ),
                        name="3σ position ellipsoid",
                    )
                ],
                traces=[6],
            )
        )

    max_extent = max(
        1.0,
        max(
            float(np.max(np.abs(points)))
            for frame in frames
            for points in (frame.data[0].x, frame.data[0].y, frame.data[0].z)
        ),
    )
    fig.add_trace(
        go.Scatter3d(
            x=[0.0],
            y=[0.0],
            z=[0.0],
            mode="markers",
            marker={"color": AMBER, "size": 4},
            name="Estimate center",
            hovertemplate="Estimate center<br>Local ECI offset = 0<extra></extra>",
        ),
        row=3,
        col=1,
    )
    fig.add_trace(
        go.Surface(
            x=initial_surface[0],
            y=initial_surface[1],
            z=initial_surface[2],
            colorscale=[[0.0, CYAN], [1.0, ORANGE]],
            opacity=0.72,
            showscale=False,
            name="3σ position ellipsoid",
        ),
        row=3,
        col=1,
    )
    # The animated surface trace is index 7 after six variance traces and the center marker.
    frames = [
        go.Frame(name=frame.name, data=frame.data, traces=[7])
        for frame in frames
    ]
    fig.frames = frames

    slider_steps = [
        {
            "args": [
                [frame.name],
                {
                    "mode": "immediate",
                    "frame": {"duration": 0, "redraw": True},
                    "transition": {"duration": 0},
                },
            ],
            "label": f"{float(frame.name):.1f}",
            "method": "animate",
        }
        for frame in frames
    ]
    fig.update_layout(
        title={
            "text": "ORBITAL TELEMETRY  //  NAVIGATION COVARIANCE",
            "x": 0.035,
            "font": {"size": 22, "color": PEARL, "family": "monospace"},
        },
        template="plotly_dark",
        paper_bgcolor=BACKGROUND,
        plot_bgcolor=PANEL,
        font={"color": TEXT, "family": "monospace"},
        height=1320,
        margin={"l": 70, "r": 40, "t": 95, "b": 80},
        hoverlabel={"bgcolor": PANEL, "bordercolor": PEARL, "font": {"color": TEXT}},
        scene={
            "bgcolor": BACKGROUND,
            "xaxis": {
                "title": "Local ΔX (m)",
                "range": [-max_extent, max_extent],
                "backgroundcolor": PANEL,
                "gridcolor": GRID,
            },
            "yaxis": {
                "title": "Local ΔY (m)",
                "range": [-max_extent, max_extent],
                "backgroundcolor": PANEL,
                "gridcolor": GRID,
            },
            "zaxis": {
                "title": "Local ΔZ (m)",
                "range": [-max_extent, max_extent],
                "backgroundcolor": PANEL,
                "gridcolor": GRID,
            },
            "aspectmode": "data",
        },
        updatemenus=[
            {
                "type": "buttons",
                "showactive": False,
                "x": 0.02,
                "y": 0.02,
                "buttons": [
                    {
                        "label": "▶ Play",
                        "method": "animate",
                        "args": [
                            None,
                            {
                                "fromcurrent": True,
                                "frame": {"duration": 100, "redraw": True},
                                "transition": {"duration": 0},
                            },
                        ],
                    },
                    {
                        "label": "⏸ Pause",
                        "method": "animate",
                        "args": [
                            [None],
                            {
                                "mode": "immediate",
                                "frame": {"duration": 0, "redraw": False},
                                "transition": {"duration": 0},
                            },
                        ],
                    },
                ],
            }
        ],
        sliders=[
            {
                "active": 0,
                "x": 0.22,
                "y": 0.02,
                "len": 0.75,
                "currentvalue": {"prefix": "Time (s): ", "font": {"color": AMBER}},
                "steps": slider_steps,
            }
        ],
    )
    fig.update_annotations(font={"color": ORANGE, "size": 12, "family": "monospace"})
    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Generate interactive orbit, state, and covariance plots."
    )
    parser.add_argument("--input", type=Path, default=Path("output.csv"))
    parser.add_argument("--output", type=Path, default=Path("output.html"))
    args = parser.parse_args()

    required_columns = {
        "t_sec",
        "r_x_m",
        "r_y_m",
        "r_z_m",
        "v_x_m_s",
        "v_y_m_s",
        "v_z_m_s",
        "estimate_r_x_m",
        "estimate_r_y_m",
        "estimate_r_z_m",
        *covariance_columns(),
    }
    df = pd.read_csv(args.input)
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(
            f"{args.input} is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )
    if df.empty:
        raise ValueError(f"{args.input} contains no simulation samples")

    state_fig = build_state_figure(df)
    covariance_fig = build_covariance_figure(df)
    state_html = state_fig.to_html(
        full_html=False,
        include_plotlyjs=True,
        config=PLOTLY_CONFIG,
        div_id="state-plot",
    )
    covariance_html = covariance_fig.to_html(
        full_html=False,
        include_plotlyjs=False,
        config=PLOTLY_CONFIG,
        div_id="covariance-plot",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>AOCS Orbital Telemetry</title>
<style>
body {{ margin: 0; background: {BACKGROUND}; color: {TEXT}; font-family: monospace; }}
.tabs {{ display: flex; gap: 8px; padding: 12px 18px; background: {PANEL}; }}
.tab-button {{ border: 1px solid {GRID}; background: {BACKGROUND}; color: {TEXT};
  padding: 10px 16px; font: inherit; cursor: pointer; }}
.tab-button.active {{ border-color: {ORANGE}; color: {AMBER}; }}
.tab-content {{ display: none; }}
.tab-content.active {{ display: block; }}
</style>
</head>
<body>
<nav class="tabs">
  <button class="tab-button active" onclick="showTab('state-tab', this)">STATE &amp; ORBIT</button>
  <button class="tab-button" onclick="showTab('covariance-tab', this)">COVARIANCE</button>
</nav>
<main>
  <section id="state-tab" class="tab-content active">{state_html}</section>
  <section id="covariance-tab" class="tab-content">{covariance_html}</section>
</main>
<script>
function showTab(tabId, button) {{
  document.querySelectorAll('.tab-content').forEach((tab) => tab.classList.remove('active'));
  document.querySelectorAll('.tab-button').forEach((tabButton) => tabButton.classList.remove('active'));
  document.getElementById(tabId).classList.add('active');
  button.classList.add('active');
  window.setTimeout(() => {{
    document.querySelectorAll('#' + tabId + ' .js-plotly-plot').forEach((plot) => {{
      Plotly.Plots.resize(plot);
    }});
  }}, 50);
}}
</script>
</body>
</html>
""",
        encoding="utf-8",
    )
    print(f"Wrote {args.output} with {len(df)} AOCS samples")


if __name__ == "__main__":
    main()
