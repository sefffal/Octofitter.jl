
# Calculate a list of nicely spaced date ticks given a vector
# of MJD.
# Returns the locations of the ticks in MJD and the formatted strings.
function _date_ticks(ts)

    # For secondary date axis on top
    date_start_day = mjd2date(ts[begin])
    # Remove any dates before 0 BC, the Julia datetime parser
    # fails on negative years...
    date_start = max.(Date("0000-01-01"), date_start_day)

    date_end_day = mjd2date(ts[end])
    date_start = Date(Dates.year(date_start), month(date_start))

    date_end = Date(Dates.year(date_end_day), month(date_end_day))
    dates = range(date_start, date_end, step=Year(1))
    dates_str = string.(year.(dates))
    if length(dates) == 1
        date_start_month_floor = Date(Dates.year(date_start), Dates.month(date_start))
        dates = range(date_start_month_floor, date_end, step=Month(1))
        dates_str = map(d->string(Dates.year(d),"-",lpad(month(d),2,'0')),dates)
        minor_ticks = range(Date(date_start_day)+Day(1), Date(date_end_day)-Day(1), step=Day(1))
    else
        year_step = 1
        date_start_year_floor = Date(Dates.year(date_start))
        while length(dates) > 8
            year_step += 1
            dates = range(date_start_year_floor, date_end, step=Year(year_step))
        end
        dates_str = string.(year.(dates))
        if step(dates) == Year(1)
            minor_ticks = range(date_start, date_end, step=Month(1))
        else
            minor_ticks = Float64[]
        end
    end


    return (mjd.(string.(dates)), dates_str, mjd.(string.(minor_ticks)))
end


# From PairPlots.jl
function margin_confidence_default_formatter(low,mid,high)
    largest_error = max(abs(high), abs(low))
    # Fallback for series with no variance
    if largest_error == 0
        if mid == 0
            digits_after_dot = 0
        else
            digits_after_dot = max(0, 1 - round(Int, log10(abs(mid))))
        end
        @static if VERSION >= v"1.10"
            title = @sprintf(
                "%.*f",
                digits_after_dot, mid,
            )
        else
            title = @eval @sprintf(
                $("\$%.$(digits_after_dot)f\$"),
                $mid,
            )
        end
        return title
    end

    digits_after_dot = max(0, 1 - round(Int, log10(largest_error)))
    use_scientific = digits_after_dot > 4

    if use_scientific
        if round(low, digits=digits_after_dot) == round(high, digits=digits_after_dot)
            title = @sprintf(
                "\$(%.1f \\pm %.1f)\\times 10^{-%d}\$",
                mid*10^(digits_after_dot-1),
                high*10^(digits_after_dot-1),
                (digits_after_dot-1)
            )
        else
            title = @sprintf(
                "\$(%.1f^{+%.1f}_{-%.1f})\\times 10^{-%d}\$",
                mid*10^(digits_after_dot-1),
                high*10^(digits_after_dot-1),
                low*10^(digits_after_dot-1),
                (digits_after_dot-1)
            )
        end
    else
        # '*' format specifier only supported in Julia 1.10+
        @static if VERSION >= v"1.10"
            if round(low, digits=digits_after_dot) == round(high, digits=digits_after_dot)
                title = @sprintf(
                    "\$%.*f \\pm %.*f\$",
                    digits_after_dot, mid,
                    digits_after_dot, high,
                )
            else
                title = @sprintf(
                    "\$%.*f^{+%.*f}_{-%.*f}\$",
                    digits_after_dot, mid,
                    digits_after_dot, high,
                    digits_after_dot, low
                )
            end
        else
            if round(low, digits=digits_after_dot) == round(high, digits=digits_after_dot)
                title = @eval @sprintf(
                    $("\$%.$(digits_after_dot)f \\pm %.$(digits_after_dot)f\$"),
                    $mid,
                    $high,
                )
            else
                title = @eval @sprintf(
                    $("\$%.$(digits_after_dot)f^{+%.$(digits_after_dot)f}_{-%.$(digits_after_dot)f}\$"),
                    $mid,
                    $high,
                    $low
                )
            end
        end
    end

    return title 
end

concat_with_nan(mat) =
    reduce((row_A, row_B) -> [row_A; NaN; row_B], eachrow(mat), init=Float64[])


# Helper from AlgebraOfGraphics -- prevent computing layout updates after adding
# each series. This speeds up figure creation and can prevent StackOverflow errors
# from circular Observable dependencies during axis limit adjustments.
get_layout(gl::Makie.GridLayout) = gl
get_layout(f::Union{Makie.Figure, Makie.GridPosition}) = f.layout
get_layout(l::Union{Makie.Block, Makie.GridSubposition}) = get_layout(l.parent)

function update(f, fig)
    layout = get_layout(fig)
    block_updates = layout.block_updates
    layout.block_updates = true
    output = f(fig)
    layout.block_updates = block_updates
    block_updates || Makie.GridLayoutBase.update!(layout)
    return output
end
