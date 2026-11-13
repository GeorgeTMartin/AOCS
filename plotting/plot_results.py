import pandas as pd
import plotly.express as px

df = pd.read_csv("output.csv")
fig = px.scatter_3d(df, x="r_x_m", y="r_y_m", z="r_z_m")
fig.write_html("output.html")
