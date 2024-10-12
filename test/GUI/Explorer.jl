using JuChrom
using Test

import GLMakie: Figure, GLFW, RGBA
import GLMakie: closeall


############################################################################################
# JuChrom.Explorer.explorer()
############################################################################################
@testset "JuChrom.Explorer.ExplorerData()" begin
    closeall()  # close any windows that may be open

    # Default initialization
    e = JuChrom.Explorer.explorer()
    @test isa(e.fig, Figure)  # check if the figure is created
    @test isa(e.window, GLFW.Window)  # check if the window is created
    closeall()

    # Custom Figure Size
    custom_size = (800, 600)
    e = JuChrom.Explorer.explorer(figure_size=custom_size)
    @test e.fig.scene.viewport.val.widths == [custom_size...]
    closeall()

    # Custom reader function
    # function dummy_reader() 
    #     return "data" 
    # end
    # explorer = ExplorerData(reader=dummy_reader)
    # @test explorer.reader !== nothing  # Check if the reader is assigned.

    # Default background color
    e = JuChrom.Explorer.explorer()
    @test e.fig.scene.backgroundcolor.val == RGBA{Float32}(1, 1, 1)
    closeall()

    # Focus on show flag
    # Can't directly test focus, but make sure the flag doesn't prevent window creation.
    e = JuChrom.Explorer.explorer(focus_on_show=false)
    @test isa(e.window, GLFW.Window)  # window should still be created
    closeall()

    # Invalid figure size inputs
    @test_throws ArgumentError JuChrom.Explorer.explorer(figure_size=(0, 1)) 
    @test_throws ArgumentError JuChrom.Explorer.explorer(figure_size=(1, 0)) 
    closeall()
end
