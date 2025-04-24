using JuChrom
using Test

import GLMakie

############################################################################################
# JuChrom.Explorer.ExplorerData()
############################################################################################
@testset "JuChrom.Explorer.ExplorerData()" begin
    GLMakie.closeall()  # close any windows that may be open

    # Default initialization
    e = JuChrom.Explorer.ExplorerData()
    @test isa(e.fig, GLMakie.Figure)  # check if the figure is created
    @test isa(e.window, GLMakie.GLFW.Window)  # check if the window is created
    GLMakie.closeall()

    # Custom Figure Size
    custom_size = (800, 600)
    e = JuChrom.Explorer.ExplorerData(figure_size=custom_size)
    @test e.fig.scene.viewport.val.widths == [custom_size...]
    GLMakie.closeall()

    # Custom reader function
    # function dummy_reader() 
    #     return "data" 
    # end
    # explorer = ExplorerData(reader=dummy_reader)
    # @test explorer.reader !== nothing  # Check if the reader is assigned.

    # Default background color
    e = JuChrom.Explorer.ExplorerData()
    @test e.fig.scene.backgroundcolor.val == GLMakie.RGBA{Float32}(1, 1, 1)
    GLMakie.closeall()

    # Focus on show flag
    # Can't directly test focus, but make sure the flag doesn't prevent window creation.
    
    e = JuChrom.Explorer.ExplorerData(focus_on_show=false)
    @test isa(e.window, GLMakie.GLFW.Window)  # window should still be created
    GLMakie.closeall()

    # Invalid figure size inputs
    GLMakie.closeall()
    @test_throws ArgumentError JuChrom.Explorer.ExplorerData(figure_size=(0, 1)) 
    @test_throws ArgumentError JuChrom.Explorer.ExplorerData(figure_size=(1, 0)) 

    GLMakie.closeall()  # close window
end


############################################################################################
# JuChrom.Explorer.fig()
############################################################################################
@testset "JuChrom.Explorer.fig()" begin
    GLMakie.closeall()
    e = JuChrom.Explorer.ExplorerData()
    @test isa(JuChrom.Explorer.fig(e), GLMakie.Figure)
    GLMakie.closeall() 
end


############################################################################################
# JuChrom.Explorer.figure_backgroundcolor!()
############################################################################################
@testset "JuChrom.Explorer.figure_backgroundcolor!()" begin
    GLMakie.closeall()
    e = JuChrom.Explorer.ExplorerData()
    JuChrom.Explorer.figure_backgroundcolor!(e, GLMakie.RGBA{Float32}(0.0f0, 0.0f0, 0.0f0))
    @test JuChrom.Explorer.figure_backgroundcolor(e).val == GLMakie.RGBA{Float32}(
        0.0f0, 0.0f0, 0.0f0)
        GLMakie.closeall() 
end


############################################################################################
# JuChrom.Explorer.figure_backgroundcolor()
############################################################################################
@testset "JuChrom.Explorer.figure_backgroundcolor()" begin
    GLMakie.closeall()
    e = JuChrom.Explorer.ExplorerData()
    @test JuChrom.Explorer.figure_backgroundcolor(e).val == GLMakie.RGBA{Float32}(
        1.0f0, 1.0f0, 1.0f0)
    GLMakie.closeall() 
end


############################################################################################
# JuChrom.Explorer.window()
############################################################################################
@testset "JuChrom.Explorer.window()" begin
    GLMakie.closeall() 
    e = JuChrom.Explorer.ExplorerData()
    @test isa(JuChrom.Explorer.window(e), GLMakie.GLFW.Window)
    GLMakie.closeall() 
end
