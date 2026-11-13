# EGM2008

EGM2008 is a fully normalized, tide-free spherical harmonic model of earth's gravitational field.

Parameters taken from the original `.gfc` file.  
```
earth_gravity_constant      0.3986004415E+15  
radius                      0.63781363E+07  
max_degree                  2190  
errors                      calibrated  
norm                        fully_normalized  
tide_system                 tide_free  
```

Literature Source: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2011JB008916  
Data Source: https://icgem.gfz.de/tom_longtime  

The original `EGM2008.gfc` file is retained at the repository root. Its checked-in binary, `data/earth_gravity/emg2008_coefficients.bin`, stores the model metadata and all 2,401,334 coefficient records. Each little-endian record contains degree (row) and order (column) as 16-bit unsigned integers, followed by C, S, sigma C, and sigma S as 64-bit doubles. The binary is tracked with Git LFS.

The simulation can load and query the dataset through `EGM2008Coefficients`:

```cpp
EGM2008Coefficients model("data/earth_gravity/emg2008_coefficients.bin");
const auto& c20 = model.coefficient(2, 0);
```

The simulation evaluates the fully normalized spherical-harmonic expansion using the tide-free EGM2008 coefficients and Vallado's spherical-coordinate acceleration formulation. It uses the static field only; no dynamic solid-Earth tide corrections are applied. The point-mass model remains available as an alternative.

Select the model in `sim/config.json` with `environment.gravity.model` (`"egm2008"` or `"point_mass"`). For EGM2008, set `environment.gravity.degree` and `environment.gravity.order` to non-negative integer truncation limits; order cannot exceed degree, and degree cannot exceed the dataset maximum. These limits affect the spherical-harmonic sum. The field remains static and tide-free, with no dynamic solid-Earth tide corrections.