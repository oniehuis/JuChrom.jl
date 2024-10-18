using JuChrom
using Test

import GLMakie
import GLMakie: Axis, Point, Screen
import GLMakie: closeall, display, events, ispressed, mouseposition
import Makie: Keyboard, KeyEvent, Mouse, MouseButtonEvent
import Makie: colorbuffer, rich, theme


############################################################################################
# JuChrom.axisbounds(axis::Axis)
############################################################################################
@testset "axisbounds(axis::Axis)" begin
    closeall()  # close any windows that may be open

    fig = GLMakie.Figure(size=(100, 100))
    screen = display(Screen(visible = false), fig)
    colorbuffer(screen)
    ax = Axis(fig[1, 1])
    
    @test JuChrom.axisbounds(ax) == (41, 84, 39, 84)
    @test_throws MethodError JuChrom.axisbounds(fig)

    closeall()  # close window
end


############################################################################################
# JuChrom.decompose_exponential(x::Float64)
############################################################################################
@testset "decompose_exponential Tests" begin
    # Test for a positive number
    @testset "Positive Number" begin
        mantissa, exponent = JuChrom.decompose_exponential(2.323 * 10^4)
        @test mantissa ≈ 2.323
        @test exponent == 4
    end

    # Test for a negative number
    @testset "Negative Number" begin
        mantissa, exponent = JuChrom.decompose_exponential(-2.323 * 10^4)
        @test mantissa ≈ -2.323
        @test exponent == 4
    end

    # Test for a small positive number
    @testset "Small Positive Number" begin
        mantissa, exponent = JuChrom.decompose_exponential(0.00123)
        @test mantissa ≈ 1.23
        @test exponent == -3
    end

    # Test for a small negative number
    @testset "Small Negative Number" begin
        mantissa, exponent = JuChrom.decompose_exponential(-0.00123)
        @test mantissa ≈ -1.23
        @test exponent == -3
    end

    # Test for zero
    @testset "Zero" begin
        mantissa, exponent = JuChrom.decompose_exponential(0.0)
        @test mantissa ≈ 0.0
        @test exponent == 0
    end

    # Test for a large positive number
    @testset "Large Positive Number" begin
        mantissa, exponent = JuChrom.decompose_exponential(1.0e10)
        @test mantissa ≈ 1.0
        @test exponent == 10
    end
    
    # Test for a large negative number
    @testset "Large Negative Number" begin
        mantissa, exponent = JuChrom.decompose_exponential(-1.0e10)
        @test mantissa ≈ -1.0
        @test exponent == 10
    end
end


############################################################################################
# JuChrom.height(ax::Axis) 
############################################################################################
@testset "JuChrom.height(ax::Axis)" begin
    closeall()  # close any windows that may be open

    fig = GLMakie.Figure(size=(0, 100))
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

    fig = GLMakie.Figure(size=(100, 100))
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

    fig = GLMakie.Figure(size=(100, 100))
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

    fig = GLMakie.Figure(size=(100, 100))
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

    fig = GLMakie.Figure(size=(100, 100))
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

    fig = GLMakie.Figure(size=(100, 100))
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

    fig = GLMakie.Figure(size=(100, 100))
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

    fig = GLMakie.Figure(size=(100, 100))
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
# JuChrom.textdimensions(text::AbstractString; font=theme(:fonts)[:regular].val, 
# fontsize=theme(:fontsize).val)
############################################################################################
@testset "textdimensions Tests" begin
    # Test for a simple text with standard font and fontsize
    @testset "Standard Font and Fontsize" begin
        width, height = JuChrom.textdimensions("Hello, world!", 
            font="TeX Gyre Heros Makie", fontsize=14)
        @test width ≈ 77.01399397850037
        @test height ≈ 16.309999465942383
    end

    # Test for a simple text with a larger fontsize
    @testset "Larger Fontsize" begin
        width, height = JuChrom.textdimensions("Hello, world!", 
            font="TeX Gyre Heros Makie", fontsize=28)
        @test width ≈ 154.02798795700073
        @test height ≈ 32.619998931884766
    end
    
    # Test for a simple text with a smaller Float64 fontsize
    @testset "Smaller Fontsize" begin
        width, height = JuChrom.textdimensions("Hello, world!", 
            font="TeX Gyre Heros Makie", fontsize=7.5)
        @test width ≈ 41.25750058889389
        @test height ≈ 8.737499713897705
    end
    
    # Test for a different font
    @testset "Different Font" begin
        width, height = JuChrom.textdimensions("Hello, world!", 
            font="TeX Gyre Heros Makie Bold Italic", fontsize=14)
        @test width ≈ 84.01399743556976
        @test height ≈ 16.309999465942383
    end
    
    # Test for an empty string
    @testset "Empty String" begin
        width, height = JuChrom.textdimensions("", font="TeX Gyre Heros Makie", 
            fontsize=14)
        @test width ≈ 0.0
        @test height ≈ 0.0
    end
    
    # Test for a multiline text
    @testset "Multiline Text" begin
        width, height = JuChrom.textdimensions("Hello,\nworld!", 
            font="TeX Gyre Heros Makie", fontsize=14)
        @test width ≈ 37.33799910545349
        @test height ≈ 32.619998931884766
    end
end


############################################################################################
# JuChrom.textdimensions(text::Makie.RichText; font=theme(:fonts)[:regular].val, 
# fontsize=theme(:fontsize).val)
############################################################################################
@testset "textdimensions Tests" begin
    # Test for a simple text with standard font and fontsize
    @testset "Standard Font and Fontsize" begin
        width, height = JuChrom.textdimensions((rich("Hello, world!")), 
            font="TeX Gyre Heros Makie", fontsize=14)
        @test width ≈ 76.51399397850037
        @test height ≈ 15.809999346733093
    end

    # Test for a simple text with a larger fontsize
    @testset "Larger Fontsize" begin
        width, height = JuChrom.textdimensions((rich("Hello, world!")), 
            font="TeX Gyre Heros Makie", fontsize=28)
        @test width ≈ 153.52798795700073
        @test height ≈ 32.11999869346619
    end
    
    # Test for a simple text with a smaller Float64 fontsize
    @testset "Smaller Fontsize" begin
        width, height = JuChrom.textdimensions((rich("Hello, world!")), 
            font="TeX Gyre Heros Makie", fontsize=7.5)
        @test width ≈ 40.75750058889389
        @test height ≈ 8.237499684095383
    end
    
    # Test with a different font
    @testset "Different Font" begin
        width, height = JuChrom.textdimensions(rich("Hello, world!"), 
            font="TeX Gyre Heros Makie Bold Italic", fontsize=14)
        @test width ≈ 83.51399743556976
        @test height ≈ 15.809999346733093
    end
    
    # Test for an empty string
    @testset "Empty String" begin
        width, height = JuChrom.textdimensions(rich(""), font="TeX Gyre Heros Makie", 
            fontsize=14)
        @test width ≈ 0.0
        @test height ≈ 0.0
    end
    
    # Test for a multiline text
    @testset "Multiline Text" begin
        width, height = JuChrom.textdimensions(rich("Hello,\nworld!"), 
            font="TeX Gyre Heros Makie", fontsize=14)
        @test width ≈ 36.83799910545349
        @test height ≈ 35.80999934673309
    end
end


############################################################################################
# JuChrom.ticks(start::Real, stop::Real; count::Integer=5)
############################################################################################
@testset "ticks Tests" begin
    # Test for a simple range with default count
    @test JuChrom.ticks(0.0, 10.0) == 0.0:2.0:10.0

    # Test for a simple range with specified count
    @test JuChrom.ticks(0.0, 10.0, count=10) == 0.0:1.0:10.0
    
    # Test for a simple range with another specified count
    @test JuChrom.ticks(0.0, 10.0, count=3) == 0.0:5.0:10.0
    
    # Test for a range with negative values
    @test JuChrom.ticks(-10.0, 10.0) == -10.0:5.0:10.0
    
    # Test for a small range
    @test JuChrom.ticks(0.0, 1.0) == 0.0:0.2:1.0
    
    # Test for a large range
    @test JuChrom.ticks(0.0, 1000.0) == 0.0:200.0:1000.0

    # Test for a large range
    @test JuChrom.ticks(0.0, 0.0001, 
        count=13) ≈ 0.0:9.999999999999999e-6:9.999999999999999e-5
    
    # Test for a range with non-zero start
    @test JuChrom.ticks(5.0, 15.0) == 6.0:2.0:14.0
    
    # Test for a range with non-zero start and end
    @test JuChrom.ticks(3.2, 17.8) == 5.0:5.0:15.0
    
    # Test for a range with a single tick
    @test_throws ArgumentError JuChrom.ticks(0.0, 0.0)
    
    # Test for a range with a single tick and non-zero start
    @test_throws ArgumentError JuChrom.ticks(5.0, 5.0)
end


############################################################################################
# JuChrom.width(ax::Axis) 
############################################################################################
@testset "JuChrom.width(ax::Axis)" begin
    closeall()  # close any windows that may be open

    fig = GLMakie.Figure(size=(0, 100))
    screen = display(Screen(visible = false), fig)
    colorbuffer(screen)
    ax = Axis(fig[1, 1])
    
    x_start, x_stop, _, _ = JuChrom.scenebounds(ax.scene)
    @test JuChrom.width(ax).val == x_stop - x_start
    @test_throws MethodError JuChrom.width(fig)
    
    closeall()  # close window
end
