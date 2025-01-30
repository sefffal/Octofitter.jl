# Orbit Visualization with `octoplot`
octoplot is a versatile visualization function that creates publication-quality figures showing orbit fits and data. It can generate multi-panel figures combining:

* Projected orbits in the plane of the sky (astrometry)
* Physical orbits in AU
* Time series of separations and position angles
* Radial velocity curves
* Proper motion anomaly
* Mass posteriors
* And more

Here's a basic example showing how to create a plot from your MCMC chain:
```julia
using Octofitter
# After running your MCMC fit...
fig = octoplot(model, chain)
```

By default, `octoplot` will automatically detect what kinds of data are present in your model and create appropriate panels. For example, if you have both astrometry and radial velocity data, it will show both an orbit plot and an RV curve.
You can control which panels appear using boolean flags:
```julia
fig = octoplot(model, chain;
    show_astrom=true,         # Show orbit in sky plane (mas)
    show_physical_orbit=true, # Show orbit in physical units (AU)
    show_astrom_time=true,    # Show sep/PA vs time
    show_rv=true,             # Show stellar RV
    show_relative_rv=true,    # Show planet-star relative RV
    show_pma=true,            # Show proper motion anomaly
    show_mass=true            # Show mass posterior
)
```

By default, `octoplot` draws 250 orbits randomly from your posterior samples. You can adjust this using the `N` parameter:
```julia
# Plot fewer orbits for faster rendering
fig = octoplot(model, chain, N=50)

# Plot specific posterior samples
fig = octoplot(model, chain, ii=[1,2,3])  # Plot first three samples
idx_MAP = argmax(chain[:logpost])
fig = octoplot(model, chain, ii=[idx_MAP])  # Plot maximum posterior sample
```


# Panel Types

`octoplot` creates a vertical stack of panels based on the data present in your model and the display options you select. The panels share consistent formatting - all time-based plots share aligned time axes, orbits for each planet use consistent colors, and epoch markers (if specified) appear consistently across all applicable panels.

## Astrometry Panels

### Sky-Projected Orbits (`show_astrom=true`)
Shows orbits projected onto the plane of the sky, with ΔRA and ΔDec measured in milliarcseconds (mas). The central star is marked with a star symbol at (0,0). If you have relative astrometry data, it will be plotted with error bars. The color of each orbit indicates the mean anomaly (orbit phase), progressing from periastron.

### Physical Orbits (`show_physical_orbit=true`) 
Similar to the sky-projected plot, but shows orbits in their true physical scale measured in astronomical units (AU). This can be helpful for understanding the true geometry of the system, especially for orbits viewed at high inclination.

## Time Series Panels

### Astrometry vs Time (`show_astrom_time=true`)
Two linked panels showing:
- Projected separation vs time (top)
- Position angle vs time (bottom)

If you specified `mark_epochs_mjd`, the predicted separation and position angle at those epochs will be marked. 

### Radial Velocity (`show_rv=true`)
Shows the stellar radial velocity curve(s). If you have data from multiple instruments, each will be plotted in its own panel. The model includes:
- Raw RV measurements with error bars
- lines showing individual orbit draws from the posterior
- Colored bands showing uncertainty including jitter and GP model (if present)

If you specified `mark_epochs_mjd`, the predicted RV at those epochs will be marked. See `rvpostplot` for another way to plot RV data.

### Relative Radial Velocity (`show_relative_rv=true`)
Shows the relative radial velocity between the planet and star. This panel appears when you have `PlanetRelativeRVLikelihood` data in your model.

### Proper Motion Anomaly (`show_hgca=true`)
Multiple panels showing proper motion data from the Hipparcos-Gaia Catalogue of Accelerations (HGCA):
- Proper motion in RA vs time
- Proper motion in Dec vs time
- 2D Proper motion residual plots at the Hipparcos, Gaia, and Hipparcos-Gaia epochs.

### Mass Posterior (`show_mass=true`)
A mini corner plot showing the mass posterior(s) for your planet(s). This panel appears at the bottom of the figure when mass is a parameter in your model.

Each of these panels will only appear if:
1. The relevant display option is set to `true`
2. Your model includes the appropriate type of data/likelihood
3. Your model parameterization supports that type of visualization (e.g., `show_physical_orbit` requires a full 3D orbit parameterization)


# Customizing Appearance

## Colormaps
The default colormap ("plasma") is used to indicate orbital phase in most panels, varying smoothly from periastron through the orbit. You can customize this in several ways:

```julia
# Use a different colormap from the ColorSchemes package
octoplot(model, chain, colormap=Makie.cgrad(:viridis))

# Use a gradient from light grey to a specific color
octoplot(model, chain, colormap="#0072b2")  # Blue to grey
```

For models with multiple planets, `octoplot` automatically assigns a different base color to each planet's orbit in the astrometry and astrometry-time panels. The other panels continue to use the default colormap to show orbital phase.

!!! tip
    When choosing a colormap, avoid categorical colormaps in favor of those that vary smoothly. If you don't need to mark the exact location of periastron, you can also use a cyclical colormap.


## Transparency
The `alpha` parameter controls the transparency of orbit lines. By default, it is automatically scaled based on the number of orbits being plotted to prevent overplotting:
```julia
# Override the default transparency
octoplot(model, chain, alpha=0.1)  # More transparent
```


## Figure Scale
The overall size of the figure can be adjusted using the `figscale` parameter:
```julia
# Make the figure 50% larger
octoplot(model, chain, figscale=1.5)
```

## Time Range
You can control the time span of the orbital plots using the ts parameter:
```julia
# Custom time range in Modified Julian Days
ts = range(50000, 55000, length=200)  # 200 points between MJD 50000 and 55000
octoplot(model, chain, ts=ts)
```

!!! note
    The number of points should be roughly 150 per orbital period displayed to ensure smooth curves. For highly eccentric orbits, you may need more points to capture the rapid motion near periastron.

## Marking Specific Epochs
You can highlight specific epochs across all applicable panels using mark_epochs_mjd:
```julia
# Mark three specific dates
epochs = [
    mjd("2024-01-01"),
    mjd("2025-01-01"),
    mjd("2026-01-01")
]
octoplot(model, chain, mark_epochs_mjd=epochs)
```
These markers appear consistently across all panels, using the same color and style to show model predictions at those specific times.

## Post-Creation Customization
Since octoplot returns a Makie figure object, you can further customize the plot after creation:
```julia
# Create the plot
fig = octoplot(model, chain)

# Access and modify specific axes
ax_orbit = fig.content[1]  # First axis (usually the orbit plot)
xlims!(ax_orbit, -100, 100)  # Set x-axis limits in mas
ylims!(ax_orbit, -100, 100)  # Set y-axis limits in mas

# Add a title
ax_orbit.title = "HD 12345 Orbital Fit"
fig
```

# Time Range Control

The time span shown in orbit plots can be controlled using the `ts` parameter. By default, `octoplot` automatically determines an appropriate range based on:
- Your data epochs
- The median orbital period from your posterior
- A small padding factor for visual clarity

You can override this behavior by providing your own time range:

```julia
# Custom time range
ts = range(mjd("2020-01-01"), mjd("2025-01-01"), length=200)
octoplot(model, chain, ts=ts)
```

!!! note
    The `ts` parameter only affects time-based panels (RV curves, proper motion, etc). The sky-projected orbit plots (`show_astrom`) and physical orbit plots (`show_physical_orbit`) use a separate internal algorithm to ensure smooth curves.

!!! tip
    If you have widely separated data epochs, you might want to zoom in on specific time ranges to better see the detail in your data. For example:

```julia
# Zoom in on first epoch
ts = range(mjd("2020-01-01"), mjd("2021-01-01"), length=200)
fig1 = octoplot(model, chain, ts=ts)
```


# Post-Creation Customization

Since `octoplot` returns a Makie figure object, you can further customize any aspect of the plot after creation. The figure contains a vertical layout of different plot panels depending on which elements you chose to display.

## Accessing Plot Elements
```julia
# Create the plot
fig = octoplot(model, chain)

# Access the axes
axes = fig.content    # Get all axes

# Common panel indices
orbit_ax = fig.content[1]     # Sky-projected orbit plot (if show_astrom=true)
rv_ax = fig.content[2]        # RV plot (if show_rv=true)
# etc.
```


## Common Adjustments
```julia
# Adjust axis limits
xlims!(orbit_ax, -100, 100)   # Set x-axis limits in mas
ylims!(orbit_ax, -100, 100)   # Set y-axis limits in mas

# Change axis labels
orbit_ax.xlabel = "ΔRA [mas]"
orbit_ax.ylabel = "ΔDec [mas]"

# Add a title
orbit_ax.title = "HD 12345 Orbital Fit"

# Adjust legend
Legend(fig[1,2], orbit_ax, "Posterior Draws")

# Save the modified figure
save("orbit_plot.pdf", fig)
```

!!! tip
    Take care when modifying time-based panels (RV, proper motion, etc) as they share synchronized x-axes. Modifying the time limits of one panel will affect all others.


# Saving Figures

When using CairoMakie (recommended for publication-quality outputs), you can save figures in several formats:

```julia
fig = octoplot(model, chain)

# Save as PNG (default)
save("orbit_plot.png", fig)

# Save as PDF (great for publications)
save("orbit_plot.pdf", fig)

# Save as SVG (good for further editing)
save("orbit_plot.svg", fig)

# Increase PNG resolution
save("orbit_plot.png", fig, px_per_unit=5)  # 5x default resolution
```

!!! note
    Many plot elements are internally rasterized to 4 points per pixel for performance, so extremely high px_per_unit values may not improve quality.

If using GLMakie instead of CairoMakie, you get interactive figures that can be zoomed and panned before saving. However, only PNG output is supported. GLMakie is great for exploration, while CairoMakie is preferred for final publication figures.