
# Calculate a list of nicely spaced date ticks given a vector
# of MJD.
# Returns the locations of the ticks in MJD and the formatted strings.
function _date_ticks(ts)

    # For secondary date axis on top
    date_start = mjd2date(ts[begin])
    # Remove any dates before 0 BC, the Julia datetime parser
    # fails on negative years...
    date_start = max.(Date("0000-01-01"), date_start)

    date_end = mjd2date(ts[end])
    date_start = Date(Dates.year(date_start), month(date_start))

    date_end = Date(Dates.year(date_end), month(date_end))
    dates = range(date_start, date_end, step=Year(1))
    dates_str = string.(year.(dates))
    if length(dates) == 1
        dates = range(date_start, date_end, step=Month(1))
        dates_str = map(d->string(Dates.year(d),"-",lpad(month(d),2,'0')),dates)
    else
        year_step = 1
        while length(dates) > 8
            year_step += 1
            dates = range(date_start, date_end, step=Year(year_step))
        end
        dates_str = string.(year.(dates))
    end
    return (mjd.(string.(dates)), dates_str)
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
                "\$%.*f",
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

# https://discourse.julialang.org/t/equivalent-of-matlabs-unwrap/44882/4?
function unwrap!(x, period = 2π)
	y = convert(eltype(x), period)
	v = first(x)
	@inbounds for k = eachindex(x)
		x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
	end
end
