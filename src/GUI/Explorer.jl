module Explorer

import GLMakie: Axis, Figure, GLFW, Keyboard, Mouse, Observable, RGBA, Scene
import GLMakie: events, ispressed, mouseposition, on, viewport
import GLMakie: @lift
import GLMakie

import Makie: KeyEvent, MouseButtonEvent, theme
import Makie

export axisbounds
export decompose_exponential
export explorer
export height
export inaxis
export inscene
export keyboardbutton!
export mousebutton!
export mouseposition!
export pixelspace_figure
export scenebounds
export textdimensions
export ticks
export width

include("makie_utilities.jl")
include("ExplorerData.jl")


function explorer(;
    reader=nothing,
    figure_size=(1528, 750),
    focus_on_show=true)

    ######################################################################################
    # Enable interactive features
    ######################################################################################
    e = ExplorerData(;
        reader=reader,
        figure_size=figure_size,
        focus_on_show=focus_on_show
        )
end

end  # module