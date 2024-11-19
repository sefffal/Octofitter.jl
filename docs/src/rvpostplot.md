# RV Visualization with `rvpostplot`

While `octoplot` provides a broad overview of all your data and orbital fits, `rvpostplot` specializes in detailed visualization of radial velocity data. It creates a multi-panel figure showing:
- The full RV time series with model fits
- Residuals from the model
- Phase-folded curves for each planet

Two versions are available:
- `rvpostplot(model, chain)`: Shows a single posterior sample
- `rvpostplot_animated(model, chain)`: Creates an animation cycling through different posterior samples

Here is an example:

![](assets/rv-postplot-1.png)

## Basic Usage
```julia
# Plot a single sample (by default, the maximum posterior sample)
fig = rvpostplot(model, chain)

# Create an animation
fig = rvpostplot_animated(model, chain)
```


## Understanding the Plot Panels
### Time Series Panel
The top panel shows:

* RV measurements from each instrument (different colors)
* Model fits including any Gaussian Process stellar activity model
* Error bars:
    * Colored bars: Raw measurement uncertainty
    * Grey bars: Measurement + instrument jitter
    * Colored bands: uncertainty from the GP model (if used)
* Optional perspective aka. secular acceleration line for models based on `AbsoluteVisual{...}` orbit

##  Residuals Panel
Shows the difference between the data and model.
Note that in the residuals and phase-folded plots, the grey bars included the GP uncertainty too.

## Phase-Folded Panels
For each planet in your model, a phase-folded panel shows:
* Data folded at the planet's orbital period
* Other planet signals subtracted from the data
* Binned data points (red) with uncertainties
* Model curve in blue

Optional text summary showing orbital parameters and uncertainties (pass `show_summary = true`)


## Detailed Options and Customization

### Panel Selection
```julia
# Plot the maximum posterior sample with options
rvpostplot(model, chain;
    show_perspective=true,   # Show perspective acceleration line
    show_summary=true,      # Show orbital parameter summary text
)
```

### Orbit Sample Selection
By default, `rvpostplot` shows the maximum posterior sample. You can specify a different sample:
```julia
# Plot a specific sample
i_sample = 42
fig = rvpostplot(model, chain, i_sample)

# Plot the first sample
fig = rvpostplot(model, chain, 1)
```

### Animation Options
The `rvpostplot_animated` function creates an animation that cycles through different posterior samples, helping visualize the range of orbits consistent with your data.

```julia
# Basic animation with default settings
anim = rvpostplot_animated(model, chain)

# Customize animation parameters
anim = rvpostplot_animated(model, chain;
    N = 50,            # Number of frames (default)
    framerate = 4,     # Frames per second
    compression = 1,   # Video compression level
    fname = "rv-posterior.mp4"  # Output filename
)
```
!!! note
    The default of 50 frames usually provides a good balance between smooth animation and reasonable processing time.

### Advanced Animation Control
For fine-grained control over each frame, you can provide a callback function:
```julia
# Example: Customize axis limits for each frame
function adjust_frame(fig)
    ax = fig.content[1]  # Get first axis
    ylims!(ax, -100, 100)  # Set y limits
    return fig
end

anim = rvpostplot_animated(model, chain, callback=adjust_frame)
```