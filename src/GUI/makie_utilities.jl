""" 
    JuChrom.axisbounds(axis::Axis)

Return the bounds of a 2D `axis` within the figure's pixel space.

# Examples
```jldoctest
julia> using Makie;

julia> fig = Figure(size=(100, 100));

julia> ax = Axis(fig[1, 1]);

julia> x_start, x_stop, y_start, y_stop = JuChrom.axisbounds(ax)
(41, 84, 39, 84)
```
"""
axisbounds(axis::Axis) = scenebounds(axis.scene)


"""
    JuChrom.decompose_exponential(x::Float64) -> Tuple{Float64, Int}

Decomposes a number written in scientific notation into its mantissa (the base value) and 
exponent, where the base of the exponent is 10.

# Examples
```jldoctest
julia> JuChrom.decompose_exponential(2.323 * 10^4) .≈ (2.323, 4)
(true, true)

julia> JuChrom.decompose_exponential(12345.0) .≈ (1.2345, 4)
(true, true)
```
"""
function decompose_exponential(x::Real)
    x == 0 && return 0.0, 0
    exponent::Int = floor(log10(abs(x)))
    mantissa = x / 10.0^exponent
    mantissa, exponent
end


""" 
    JuChrom.height(ax::Axis) -> Observable{Int}

Return the axis height in pixels as an observable.

# Examples
```jldoctest
julia> using Makie;

julia> fig = Figure(size=(100, 100));

julia> ax = Axis(fig[1, 1]);

julia> JuChrom.height(ax)
Observable(45)
```
"""
height(ax::Axis)::Observable{Int} = @lift($(ax.yaxis.attributes.endpoints)[2][2] 
    .- $(ax.yaxis.attributes.endpoints)[1][2])


""" 
    JuChrom.inaxis(axis::Axis, xy)

Checks whether the 2D coordinates `xy` in the figure's the pixel space are within the 
boundaries of the 2D `axis`. The function returns `true` if the point lies within the 
boundaries and `false` otherwise.

# Examples
```jldoctest
julia> using Makie;

julia> fig = Figure(size=(100, 100));

julia> ax = Axis(fig[1, 1]);

julia> x_start, x_stop, y_start, y_stop = JuChrom.scenebounds(ax.scene);

julia> JuChrom.inaxis(ax, (x_start, y_start))
true

julia> JuChrom.inaxis(ax, (x_stop, y_stop))
true

julia> JuChrom.inaxis(ax, (x_start, y_stop))
true

julia> JuChrom.inaxis(ax, (x_stop, y_start))
true

julia> JuChrom.inaxis(ax, (x_start - 1, y_start))
false

julia> JuChrom.inaxis(ax, (x_start, y_start - 1))
false
```
"""
inaxis(axis::Axis, xy) = inscene(axis.scene, xy)


""" 
    JuChrom.inscene(scene::Scene, xy)

Checks whether the 2D coordinates `xy` in the figure's the pixel space are within the 
boundaries of the 2D `scene`.  The function returns `true` if the point lies within the 
boundaries and `false` otherwise.

# Examples
```jldoctest
julia> using Makie;

julia> fig = Figure(size=(100, 100));

julia> ax = Axis(fig[1, 1]);

julia> x_start, x_stop, y_start, y_stop = JuChrom.scenebounds(ax.scene);

julia> JuChrom.inscene(ax.scene, (x_start, y_start))
true

julia> JuChrom.inscene(ax.scene, (x_stop, y_stop))
true

julia> JuChrom.inscene(ax.scene, (x_start, y_stop))
true

julia> JuChrom.inscene(ax.scene, (x_stop, y_start))
true

julia> JuChrom.inscene(ax.scene, (x_start - 1, y_start))
false

julia> JuChrom.inscene(ax.scene, (x_start, y_start - 1))
false
```
"""
function inscene(scene::Scene, xy)
    length(xy) == 2 || throw(
        ArgumentError("the mouse position is not a vector with exactly two elements"))
    x_start, x_stop, y_start, y_stop = scenebounds(scene)
    (x_start ≤ first(xy) ≤ x_stop) && (y_start ≤ last(xy) ≤ y_stop)
end


""" 
    JuChrom.keyboardbutton!(fig::Figure, event::Makie.KeyEvent)

Assign the KeyEvent to the figure's observable keyboardbutton, which maintains the 
current state of the keyboard buttons.

# Examples
```jldoctest
julia> using Makie;

julia> fig = Figure(size=(100, 100));

julia> ax = Axis(fig[1, 1]);

julia> JuChrom.keyboardbutton!(fig, Makie.KeyEvent(Keyboard.b, Keyboard.press))
Makie.KeyEvent(Makie.Keyboard.b, Makie.Keyboard.press)

julia> events(fig).keyboardstate
Set{Makie.Keyboard.Button} with 1 element:
  Makie.Keyboard.b

julia> events(fig).keyboardstate == Set([Keyboard.b])
true

julia> events(fig).keyboardstate == Set([Keyboard.p])
false

julia> ispressed(events(fig), Keyboard.b)
true

julia> ispressed(events(fig), Keyboard.p)
false

julia> JuChrom.keyboardbutton!(fig, Makie.KeyEvent(Keyboard.b, Keyboard.release))
Makie.KeyEvent(Makie.Keyboard.b, Makie.Keyboard.release)

julia> ispressed(events(fig), Keyboard.b)
false
```
"""
keyboardbutton!(fig::Figure, event::KeyEvent) = events(fig).keyboardbutton[] = event


""" 
    JuChrom.mousebutton!(fig::Figure, event::Makie.MouseButtonEvent)

Assign the MouseButtonEvent to the figure's observable mousebutton, which maintains the 
current state of the mouse buttons.

# Examples
```jldoctest
julia> using Makie;

julia> fig = Figure(size=(100, 100));

julia> ax = Axis(fig[1, 1]);

julia> JuChrom.mousebutton!(fig, Makie.MouseButtonEvent(Mouse.left, Mouse.press))
Makie.MouseButtonEvent(Makie.Mouse.left, Makie.Mouse.press)

julia> events(fig).mousebutton[]
Makie.MouseButtonEvent(Makie.Mouse.left, Makie.Mouse.press)

julia> events(fig).mousebuttonstate
Set{Makie.Mouse.Button} with 1 element:
  Makie.Mouse.left

julia> events(fig).mousebuttonstate == Set([Mouse.left])
true

julia> events(fig).mousebuttonstate == Set([Mouse.right])
false

julia> ispressed(events(fig), Mouse.left)
true

julia> ispressed(events(fig), Mouse.right)
false

julia> JuChrom.mousebutton!(fig, Makie.MouseButtonEvent(Mouse.left, Mouse.release))
Makie.MouseButtonEvent(Makie.Mouse.left, Makie.Mouse.release)

julia> ispressed(events(fig), Mouse.left)
false
```
"""
mousebutton!(fig::Figure, event::MouseButtonEvent) = events(fig).mousebutton[] = event


""" 
    JuChrom.mouseposition!(fig::Figure, pos)

Assign the coordinates pos to the figure's observable mouseposition, which stores the 
current position of the mouse in the figure's pixel space.

# Examples
```jldoctest
julia> using Makie;

julia> fig = Figure(size=(100, 100));

julia> ax = Axis(fig[1, 1]);

julia> JuChrom.mouseposition!(fig, (0, 0))
(0, 0)

julia> mouseposition(fig)
2-element Point{2, Float64} with indices SOneTo(2):
 0.0
 0.0

julia> JuChrom.mouseposition!(fig, (15.2, 20.1))
(15.2, 20.1)

julia> mouseposition(fig) ≈ Point{2, Float64}(15.199999809265137, 20.100000381469727)
true 

julia> JuChrom.mouseposition!(fig, JuChrom.pixelspace_figure(ax.scene, (10, 10)))
(51, 49)
```
"""
mouseposition!(fig::Figure, pos) = events(fig).mouseposition[] = pos


""" 
    JuChrom.pixelspace_figure(scene::Scene, xy) 

Return the coordinates in the figure's pixel space when given pixel coordinates relative 
to the scene. For example, if the scene's origin is at (50, 60) and the function receives 
the scene pixel coordinates (15, 10), it would return (65, 70).

# Examples
```jldoctest
julia> using Makie;

julia> fig = Figure(size=(100, 100));

julia> ax = Axis(fig[1, 1]);

julia> x_start, x_stop, y_start, y_stop = JuChrom.scenebounds(ax.scene)
(41, 84, 39, 84)

julia> JuChrom.pixelspace_figure(ax.scene, (0, 0))
(41, 39)

julia> JuChrom.pixelspace_figure(ax.scene, (15, 10))
(56, 49)

julia> JuChrom.pixelspace_figure(ax.scene, (-5, -5))
(36, 34)
```
"""
pixelspace_figure(scene::Scene, xy) = ((xy .+ viewport(scene)[].origin)...,)


""" 
    JuChrom.scenebounds(scene::Scene)

Return the bounds of a 2D `scene` within the figure's pixel space.

# Examples
```jldoctest
julia> using Makie;

julia> fig = Figure(size=(100, 100));

julia> ax = Axis(fig[1, 1]);

julia> x_start, x_stop, y_start, y_stop = JuChrom.scenebounds(ax.scene)
(41, 84, 39, 84)
```
"""
function scenebounds(scene::Scene)
    length(viewport(scene)[].origin) == 2 || throw(
        ArgumentError("the scene is not confined to exactly two dimensions"))
    x, y = viewport(scene)[].origin
    w, h = viewport(scene)[].widths
    x, x + w, y, y + h
end


"""
    JuChrom.textdimensions(text::AbstractString; font=theme(:fonts)[:regular].val, 
        fontsize=theme(:fontsize).val) -> NTuple{2, Float64}

Return the width and height of the given `text` in pixels. Th optional keyword arguments
`font` and `fontsize` allow you to specify the font and font size to use for the text. They 
default to the regular font and font size from the current theme.

# Example
```jldoctest
julia> using GLMakie;

julia> text = "Hello, world!"
"Hello, world!"

julia> (w, h) = JuChrom.textdimensions(text, font="TeX Gyre Heros Makie", fontsize=24);

julia> (w, h) .≈ (132.02399730682373, 27.959999084472656)
(true, true)
``` 
"""
function textdimensions(text::AbstractString; font=theme(:fonts)[:regular].val, 
        fontsize=theme(:fontsize).val)::NTuple{2, Float64}
    bbox = Makie.text_bb(text, font, fontsize)

    width = bbox.widths[1]
    height = bbox.widths[2]

    width, height
end


"""
    JuChrom.textdimensions(text::Makie.RichText; font=theme(:fonts)[:regular].val, 
        fontsize=theme(:fontsize).val) -> NTuple{2, Float64}

Return the width and height of the given rich `text` in pixels. Th optional keyword 
arguments `font` and `fontsize` allow you to specify the font and font size to use for the 
text. They default to the regular font and font size from the current theme.

# Example
```jldoctest
julia> using GLMakie;

julia> richtext = rich("Hello, world!", font="TeX Gyre Heros Makie", fontsize=24)
RichText: "Hello, world!"

julia> (w, h) = JuChrom.textdimensions(richtext);

julia> (w, h) .≈ (131.52399730682373, 27.459999084472656)
(true, true)
```
"""
function textdimensions(text::Makie.RichText; font=theme(:fonts)[:regular].val, 
        fontsize=theme(:fontsize).val)::NTuple{2, Float64}
    
    # Makie.text! throws an error when provided with an empty rich text
    if length(text.children) == 1 && text.children[1] == ""
        return 0.0, 0.0
    end

    scene = Makie.Scene(size=(1, 1))
    textplot = Makie.text!(scene, [0], [0], text=text, font=font, fontsize=fontsize)
    bbox = Makie.boundingbox(textplot, textplot.markerspace[])
    
    width = bbox.widths[1] - bbox.origin[1]
    height = bbox.widths[2] - bbox.origin[2]

    width, height
end


"""
    JuChrom.ticks(start::Real, stop::Real; count::Integer=5) -> StepRangeLen

Return a sequence of tick values for a specified range, ensuring that the intervals 
between ticks are rounded to reasonable values. The function takes the start and 
stop values of the range as arguments. Additionally, the optional `count` argument 
specifies the desired number of tick values. The actual number of ticks generated may 
vary to ensure that the intervals are evenly spaced and easy to read. Defaults to 5.

# Example
```jldoctest
julia> JuChrom.ticks(3.2, 17.8)
5.0:5.0:15.0

julia> JuChrom.ticks(3.2, 17.8, count=10)
4.0:2.0:16.0
```
"""
function ticks(start::Real, stop::Real; count::Integer=5)
    start == stop && throw(ArgumentError("start and stop must be different"))
    stepwidth = (stop - start) / (count - 1)
    magnitude = 10.0^floor(log10(stepwidth))
    norminterval = stepwidth / magnitude
    stepwidth = if norminterval < 1.5
        magnitude
    elseif norminterval < 3
        magnitude * 2
    elseif norminterval < 7
        magnitude * 5
    else
        magnitude * 10
    end
    start_adjusted = ceil(start / stepwidth) * stepwidth
    stop_adjusted = floor(stop / stepwidth) * stepwidth
    range(start_adjusted, step=stepwidth, stop=stop_adjusted)
end


""" 
    JuChrom.width(ax::Axis) -> Observable{Int}

Return the axis width in pixels as an observable.

# Examples
```jldoctest
julia> using Makie;

julia> fig = Figure(size=(100, 100));

julia> ax = Axis(fig[1, 1]);

julia> JuChrom.width(ax)
Observable(43)
```
"""
width(ax::Axis)::Observable{Int} = @lift($(ax.xaxis.attributes.endpoints)[2][1] 
    .- $(ax.xaxis.attributes.endpoints)[1][1])
