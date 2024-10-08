module Explorer

import GLMakie: Axis, Figure, GLFW, Keyboard, Mouse, Observable, RGBf, Scene
import GLMakie: events, ispressed, mouseposition, on, viewport
import GLMakie: @lift
import GLMakie

import Makie: KeyEvent, MouseButtonEvent

export axisbounds
export explorer
export height
export inaxis
export inscene
export keyboardbutton!
export mousebutton!
export mouseposition!
export pixelspace_figure
export scenebounds
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