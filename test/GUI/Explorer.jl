using JuChrom
using Test

import GLMakie: Axis, Figure, GLFW, Point, Screen, RGBA
import GLMakie: closeall, display, events, ispressed, mouseposition
import Makie: Keyboard, KeyEvent, Mouse, MouseButtonEvent
import Makie: colorbuffer


############################################################################################
# JuChrom.axisbounds(axis::Axis)
############################################################################################
@testset "axisbounds(axis::Axis)" begin
    closeall()  # close any windows that may be open

    fig = Figure(size=(100, 100))
    screen = display(Screen(visible = false), fig)
    colorbuffer(screen)
    ax = Axis(fig[1, 1])
    
    @test JuChrom.axisbounds(ax) == (41, 84, 39, 84)
    @test_throws MethodError JuChrom.axisbounds(fig)

    closeall()  # close window
end


############################################################################################
# JuChrom.Explorer.ExplorerData()
############################################################################################
@testset "JuChrom.Explorer.ExplorerData()" begin
    closeall()  # close any windows that may be open

    # Verify that ExplorerData can be initialized with default settings
    e = JuChrom.Explorer.ExplorerData()

    # Default initialization
    closeall()
    e = JuChrom.Explorer.ExplorerData()
    @test isa(e.fig, Figure)  # Check if the figure is created.
    @test isa(e.window, GLFW.Window)  # Check if the window is created.

    # Custom Figure Size
    closeall()
    custom_size = (800, 600)
    e = JuChrom.Explorer.ExplorerData(figure_size=custom_size)
    @test e.fig.scene.viewport.val.widths == [custom_size...]

    # Custom reader function
    # function dummy_reader() 
    #     return "data" 
    # end
    # explorer = ExplorerData(reader=dummy_reader)
    # @test explorer.reader !== nothing  # Check if the reader is assigned.

    # Default background color
    closeall()
    e = JuChrom.Explorer.ExplorerData()
    @test e.fig.scene.backgroundcolor.val == RGBA{Float32}(1.0f0,1.0f0,1.0f0)   

    # Focus on show flag
    # Can't directly test focus, but make sure the flag doesn't prevent window creation.
    closeall()
    e = JuChrom.Explorer.ExplorerData(focus_on_show=false)
    @test isa(e.window, GLFW.Window)  # Window should still be created.

    # Invalid figure size inputs
    # Negative figure size should throw an error
    closeall()
    @test_throws ArgumentError JuChrom.Explorer.ExplorerData(figure_size=(0, 1)) 
    @test_throws ArgumentError JuChrom.Explorer.ExplorerData(figure_size=(1, 0)) 

    closeall()  # close window
end


############################################################################################
# JuChrom.height(ax::Axis) 
############################################################################################
@testset "JuChrom.height(ax::Axis)" begin
    closeall()  # close any windows that may be open

    fig = Figure(size=(0, 100))
    screen = display(Screen(visible = false), fig)
    colorbuffer(screen)
    ax = Axis(fig[1, 1])

    _, _, y_start, y_stop = JuChrom.scenebounds(ax.scene)
    @test JuChrom.height(ax).val == y_stop - y_start
    @test_throws MethodError JuChrom.height(fig)
    
    closeall()  # close window
end


############################################################################################
# JuChrom.inaxis(axis::Axis, xy)
############################################################################################
@testset "JuChrom.inaxis(axis::Axis, xy)" begin
    closeall()  # close any windows that may be open

    fig = Figure(size=(100, 100))
    screen = display(Screen(visible = false), fig)
    colorbuffer(screen)
    ax = Axis(fig[1, 1])
    x_start, x_stop, y_start, y_stop = JuChrom.scenebounds(ax.scene)

    # Position inside the axis limits
    x = x_start + (x_stop - x_start) / 2
    y = y_start + (y_stop - y_start) / 2
    mp_inside = (x, y)  
    @test JuChrom.Explorer.inaxis(ax, mp_inside) == true

    # Position on the boundaries of the axis limits
    mp_boundary_bottom_left = (x_start, y_start)
    @test JuChrom.Explorer.inaxis(ax, mp_boundary_bottom_left) == true

    mp_boundary_top_right = (x_stop, y_stop)
    @test JuChrom.Explorer.inaxis(ax, mp_boundary_top_right) == true

    # Mouse position outside the axis limits (left and below)
    mp_left_outside = (x_start - 1, y_start)
    @test JuChrom.Explorer.inaxis(ax, mp_left_outside) == false

    mp_below_outside = (x_start, y_start - 1)
    @test JuChrom.Explorer.inaxis(ax, mp_below_outside) == false

    # Mouse position outside the axis limits (right and above)
    mp_right_outside = (x_stop + 1, y_start)
    @test JuChrom.Explorer.inaxis(ax, mp_right_outside) == false

    mp_above_outside = (x_start, y_stop + 1)
    @test JuChrom.Explorer.inaxis(ax, mp_above_outside) == false

    # Non-square axis limits
    ax_non_square = Axis(fig[1, 1], limits=(0, 800, 0, 100))
    x_start, x_stop, y_start, y_stop = JuChrom.scenebounds(ax.scene)
    x = x_start + (x_stop - x_start) / 2
    y = y_start + (y_stop - y_start) / 2
    mp_non_square = (x, y)  
    @test JuChrom.Explorer.inaxis(ax_non_square, mp_non_square) == true

    # Test: Invalid input (not exactly 2D mouse position)
    mp_too_few = (5, )  # Only 1 value
    @test_throws ArgumentError JuChrom.Explorer.inaxis(ax, mp_too_few)

    mp_too_many = (5, 5, 5)  # 3 values
    @test_throws ArgumentError JuChrom.Explorer.inaxis(ax, mp_too_many)

    # Only operated with axes
    @test_throws MethodError JuChrom.Explorer.inaxis(fig, mp_too_many)

    closeall()  # close window
end


############################################################################################
# JuChrom.inscene(scene::Scene, xy)
############################################################################################
@testset "JuChrom.inscene(scene::Scene, xy)" begin
    closeall()  # close any windows that may be open

    fig = Figure(size=(100, 100))
    screen = display(Screen(visible = false), fig)
    colorbuffer(screen)
    ax = Axis(fig[1, 1])
    x_start, x_stop, y_start, y_stop = JuChrom.scenebounds(ax.scene)

    # Position inside the axis limits
    x = x_start + (x_stop - x_start) / 2
    y = y_start + (y_stop - y_start) / 2
    mp_inside = (x, y)  
    @test JuChrom.Explorer.inscene(ax.scene, mp_inside) == true

    # Position on the boundaries of the axis limits
    mp_boundary_bottom_left = (x_start, y_start)
    @test JuChrom.Explorer.inscene(ax.scene, mp_boundary_bottom_left) == true

    mp_boundary_top_right = (x_stop, y_stop)
    @test JuChrom.Explorer.inscene(ax.scene, mp_boundary_top_right) == true

    # Mouse position outside the axis limits (left and below)
    mp_left_outside = (x_start - 1, y_start)
    @test JuChrom.Explorer.inscene(ax.scene, mp_left_outside) == false

    mp_below_outside = (x_start, y_start - 1)
    @test JuChrom.Explorer.inscene(ax.scene, mp_below_outside) == false

    # Mouse position outside the axis limits (right and above)
    mp_right_outside = (x_stop + 1, y_start)
    @test JuChrom.Explorer.inscene(ax.scene, mp_right_outside) == false

    mp_above_outside = (x_start, y_stop + 1)
    @test JuChrom.Explorer.inscene(ax.scene, mp_above_outside) == false

    # Test: Invalid input (not exactly 2D mouse position)
    mp_too_few = (5, )  # Only 1 value
    @test_throws ArgumentError JuChrom.Explorer.inscene(ax.scene, mp_too_few)

    mp_too_many = (5, 5, 5)  # 3 values
    @test_throws ArgumentError JuChrom.Explorer.inscene(ax.scene, mp_too_many)

    closeall()  # close window
end


############################################################################################
# JuChrom.keyboardbutton!(fig::Figure, event::Makie.KeyEvent)
############################################################################################
@testset "JuChrom.keyboardbutton!(fig::Figure, event::Makie.KeyEvent)" begin
    closeall()  # close any windows that may be open

    fig = Figure(size=(100, 100))
    screen = display(Screen(visible = false), fig)
    colorbuffer(screen)
    ax = Axis(fig[1, 1])

    JuChrom.keyboardbutton!(fig, KeyEvent(Keyboard.b, Keyboard.press))
    @test events(fig).keyboardstate == Set([Keyboard.b])
    @test events(fig).keyboardstate ≠ Set([Keyboard.p])
    @test ispressed(events(fig), Keyboard.b) === true
    @test ispressed(events(fig), Keyboard.p) === false

    JuChrom.keyboardbutton!(fig, KeyEvent(Keyboard.b, Keyboard.release))
    @test ispressed(events(fig), Keyboard.b) === false

    closeall()  # close window
end


############################################################################################
# JuChrom.mousebutton!(fig::Figure, event::Makie.MouseButtonEvent)
############################################################################################
@testset "JuChrom.mousebutton!(fig::Figure, event::Makie.MouseButtonEvent)" begin
    closeall()  # close any windows that may be open

    fig = Figure(size=(100, 100))
    screen = display(Screen(visible = false), fig)
    colorbuffer(screen)
    ax = Axis(fig[1, 1])

    JuChrom.mousebutton!(fig, MouseButtonEvent(Mouse.left, Mouse.press))
    @test events(fig).mousebutton[] == MouseButtonEvent(Mouse.left, Mouse.press)
    @test events(fig).mousebuttonstate == Set([Mouse.left])
    @test events(fig).mousebuttonstate ≠ Set([Mouse.right])
    @test ispressed(events(fig), Mouse.left) === true
    @test ispressed(events(fig), Mouse.right) === false

    JuChrom.mousebutton!(fig, MouseButtonEvent(Mouse.left, Mouse.release))
    @test ispressed(events(fig), Mouse.left) === false

    closeall()  # close window
end


############################################################################################
# JuChrom.mouseposition!(fig::Figure, pos)
############################################################################################
@testset "JuChrom.mouseposition!(fig::Figure, pos)" begin
    closeall()  # close any windows that may be open

    fig = Figure(size=(100, 100))
    screen = display(Screen(visible = false), fig)
    colorbuffer(screen)
    ax = Axis(fig[1, 1])

    JuChrom.mouseposition!(fig, (0, 0))
    @test mouseposition(fig) == Point{2, Float64}(0, 0)
    
    JuChrom.mouseposition!(fig, (15.2, 20.1))
    @test mouseposition(fig) ≈ Point{2, Float64}(15.199999809265137, 20.100000381469727)

    @test JuChrom.mouseposition!(fig, JuChrom.pixelspace_figure(ax.scene, (10, 10))
        ) == (51, 49)

    closeall()  # close window
end


############################################################################################
# JuChrom.pixelspace_figure(scene::Scene, xy) 
############################################################################################
@testset "JuChrom.pixelspace_figure(scene::Scene, xy)" begin
    closeall()  # close any windows that may be open

    fig = Figure(size=(100, 100))
    screen = display(Screen(visible = false), fig)
    colorbuffer(screen)
    ax = Axis(fig[1, 1])

    @test JuChrom.scenebounds(fig.scene) == (0, 100, 0, 100)
    @test JuChrom.scenebounds(ax.scene) == (41, 84, 39, 84)

    closeall()  # close window
end


############################################################################################
# JuChrom.scenebounds(scene::Scene)
############################################################################################
@testset "scenebounds(scene::Scene)" begin
    closeall()  # close any windows that may be open

    fig = Figure(size=(100, 100))
    screen = display(Screen(visible = false), fig)
    colorbuffer(screen)
    ax = Axis(fig[1, 1])

    @test JuChrom.scenebounds(fig.scene) == (0, 100, 0, 100)
    @test JuChrom.scenebounds(ax.scene) == (41, 84, 39, 84)
    @test_throws MethodError JuChrom.scenebounds(fig)
    @test_throws MethodError JuChrom.scenebounds(ax)

    closeall()  # close window
end


############################################################################################
# JuChrom.width(ax::Axis) 
############################################################################################
@testset "JuChrom.width(ax::Axis)" begin
    closeall()  # close any windows that may be open

    fig = Figure(size=(0, 100))
    screen = display(Screen(visible = false), fig)
    colorbuffer(screen)
    ax = Axis(fig[1, 1])
    
    x_start, x_stop, _, _ = JuChrom.scenebounds(ax.scene)
    @test JuChrom.width(ax).val == x_stop - x_start
    @test_throws MethodError JuChrom.width(fig)
    
    closeall()  # close window
end


# ############################################################################################
# # explorer()
# ############################################################################################
# @testset "explorer" begin
#     closeall()
    
#     fig = Figure(; size=(500, 500))
#     screen = display(Screen(visible = false), fig)
#     colorbuffer(screen)
#     ax = Axis(fig[1, 1])
#     e = events(fig)
#     e.mouseposition[] = (250, 250)
#     @test true === JuChrom.Explorer.someaction(fig, ax)
#     e.mouseposition[] = (0, 0)
#     @test false === JuChrom.Explorer.someaction(fig, ax)
# end


# fig = Figure()
# screen = display(GLMakie.Screen(visible = false), fig)
# GLMakie.Makie.colorbuffer(screen)
