# leftcolorbox(edgelength) = @lift(Rect(0, 0, $edgelength, $edgelength))
# rightcolorbox(axis, edgelength) = @lift(Rect($(size_in_px_x(axis))-$edgelength, 0, $edgelength, $edgelength))
# fontheight(font::Observable, fontsize::Observable) = @lift(height(text_bb("pÊ", $(font), $(fontsize))))

mutable struct ExplorerData
    ######################################################################################
    # Fields of the data structure.
    ######################################################################################

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Fields related to the window and its displayed figure in general.
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const fig::Figure
    figure_backgroundcolor::Observable
    const window::GLFW.Window

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 1–4 - ticklabel space that affects almost all axes equally
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # ticklabelspace_left::Observable{Float64}
    # ticklabelspace_right::Observable{Float64}

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 1 - Displays the runset name in the center and the currently set default colors
    #          for displaying data on the left and right sides of axis 2.
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # const ax1::Axis
    
    # ax1_height::Observable{Float64}
    # ax1_font::Observable{Any}
    # ax1_fontsize::Observable{Float64}
    # ax1_textcolor::Observable{Symbol}
    # ax1_run_name::Observable{String}
    # ax1_left_linecolor::Observable{Symbol}
    # ax1_right_linecolor::Observable{Symbol}

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 2 - Displays the GC data of the runsets. It consists of two parts:
    # Axis 2a: Intended to show the positions of the extracted mass spectra and of any
    #          selected diagnostic ions.
    # Axis 2b: Displays the GC data of the runsets by the two y-axes: left and right.
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 2a
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # const ax2a::Axis
    
    # ax2a_height::Observable{Float64}
    # ax2a_font::Observable{Any}
    # ax2a_fontsize::Observable{Float64}
    # ax2a_stored_data_pos_x_y::Observable{Vector{NTuple{2, Float64}}}
    # ax2a_stored_data_marker::Observable
    # ax2a_stored_data_markersize::Observable
    # ax2a_stored_data_color::Observable
    # ax2a_stored_data_plot
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 2b – Fields that do not refer or refer equally to the left and right y-axes.
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # ax2b_xaxis_required::Observable{Bool}
    # ax2b_xaxis_datatype::Observable{Symbol}
    # ax2b_xaxis_minintervalsize::Observable{Float64}
    # ax2b_xaxis_timeunit::Observable{Unitful.TimeUnits}
    # ax2b_xlims::Observable{NTuple{2, Float64}}
    # ax2b_xticks_font::Observable{Any}
    # ax2b_xticks_fontsize::Observable{Float64}
    # ax2b_yticks_font::Observable{Any}
    # ax2b_yticks_fontsize::Observable{Float64}
    # ax2b_pixelwidth::Observable{Int}
    # ax2b_pixelheight::Observable{Int}
    # ax2b_yaxes_extra_stretchfactor::Observable{Float64}
    # ax2b_stretchfactor_changespeed_slow::Float64
    # ax2b_stretchfactor_changespeed_fast::Float64
    # ax2b_linewidth_for_current_runset::Observable{Float64}
    # ax2b_linewidth_for_non_current_runset::Observable{Float64}
    # ax2b_crosshair_pos_x::Observable{Vector{Float64}}
    # ax2b_crosshair_pos_y::Observable{Vector{Float64}}
    # ax2b_crosshair_color::Observable{Symbol}
    # ax2b_crosshair_linewidth::Observable{Float64}
    # ax2b_dragstart::Vector{Float64}
    # ax2b_crosshair_label_font::Observable{Any}
    # ax2b_crosshair_label_fontsize::Observable{Float64}
    # ax2b_crosshair_label_positions::Observable{Vector{NTuple{2, Float64}}}
    # ax2b_crosshair_label_texts::Observable{Vector{String}}
    # ax2b_crosshair_label_alignments::Observable{Vector{NTuple{2, Symbol}}}
    # ax2b_crosshair_label_offsets::Observable{Vector{NTuple{2, Float64}}}
    # ax2b_crosshair_label_visible::Observable{Bool}
    # ax2b_gcmsdata_scanindexselect_start::Observable{Union{Nothing, Int}}
    # ax2b_gcmsdata_scanindexselect_stop::Observable{Union{Nothing, Int}}
    # ax2b_gcmsdata_scanindexselect_current_range::Observable{Union{Nothing, UnitRange{Int64}}}
    # ax2b_gcmsdata_scanindexselect_mode::Observable{Union{Nothing, Symbol}}
    # ax2b_gcmsdata_scanindexselect_add::Observable{Set{Int}}
    # ax2b_gcmsdata_scanindexselect_subtract::Observable{Set{Int}}
    # ax2b_gcmsdata_scantimes_add::Observable{Vector{Float64}}
    # ax2b_gcmsdata_scantimes_add_linewidth::Observable{Float64}
    # ax2b_gcmsdata_scantimes_add_color::Observable{Symbol}
    # ax2b_gcmsdata_scantimes_subtract::Observable{Vector{Float64}}
    # ax2b_gcmsdata_scantimes_subtract_linewidth::Observable{Float64}
    # ax2b_gcmsdata_scantimes_subtract_color::Observable{Symbol}

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 2b – Fields that refer to the left y-axis only.
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # const ax2b_left::Axis
    # ax2b_left_required::Observable{Bool}
    # ax2b_left_ylims::Observable{NTuple{2, Float64}}
    # ax2b_left_yaxis_stretchfactor::Observable{Float64}
    # ax2b_left_prevlimits::Vector{NTuple{4, Float64}}

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 2b – Fields that refer to the right y-axis only.
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # const ax2b_right::Axis
    # ax2b_right_required::Observable{Bool}
    # ax2b_right_ylims::Observable{NTuple{2, Float64}}
    # ax2b_right_yaxis_stretchfactor::Observable{Float64}
    # ax2b_right_prevlimits::Vector{NTuple{4, Float64}}

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 3 - Displays either a live MS or selected ion chromatograms, depending on
    #          the axis mode. It consists of two parts:
    # Axis 3a: Lists the m/z values of the selected ions whose chromatograms are
    #          displayed.
    # Axis 3b: Displays the chromatograms of the selected ions and, if the cursor is over
    #          axis 4b, the chromatogram of the ion under the cursor.
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 3a
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # const ax3a::Axis
    # ax3a_height::Observable{Float64}
    # ax3a_font::Observable{Any}
    # ax3a_fontsize::Observable{Float64}
    # ax3a_xlims::Observable{NTuple{2, Float64}}
    # ax3a_ylims::Observable{NTuple{2, Float64}}
    # ax3a_extracted_ion_strings::Observable{Vector{Makie.RichText}}
    # ax3a_extracted_ion_string_positions::Observable{Vector{NTuple{2, Float64}}}
    # ax3a_extracted_ion_string_plot
    # ax3a_ion_intensity_extrema_strings::Observable{Vector{String}}
    # ax3a_ion_intensity_extrema_string_positions::Observable{Vector{NTuple{2, Float64}}}


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 3b
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # const ax3b::Axis
    # ax3b_xlims::Observable{NTuple{2, Float64}}
    # ax3b_ylims::Observable{NTuple{2, Float64}}
    # ax3b_pixelwidth::Observable{Int}
    # ax3b_pixelheight::Observable{Int}
    # ax3b_mode::Observable{Symbol}
    # ax3b_ms_ions::Observable{Vector{Float64}}
    # ax3b_ms_lineheights::Observable{Vector{Float64}}
    # ax3b_ms_linecolor::Observable{Symbol}
    # ax3b_ms_label_positions::Observable{Vector{NTuple{2, Float32}}}
    # ax3b_ms_label_texts::Observable{Vector{String}}
    # ax3b_ion_chromatogram_scantimes::Observable{Vector{Float64}}
    # ax3b_ion_chromatogram_intensities::Observable{Vector{Float64}}
    # ax3b_ion_chromatogram_linewidth::Observable{Float64}
    # ax3b_ion_chromatogram_color::Observable{Symbol}
    # ax3b_crosshair_pos_x::Observable{Vector{Float64}}
    # ax3b_crosshair_color::Observable{Symbol}
    # ax3b_crosshair_linewidth::Observable{Float64}

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 4 - Displays the extracted MS. It consists of two parts:
    # Axis 4a: Displays the ion in the MS shown on axis 4b.
    # Axis 4b: Displays the extracted MS.
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 4a
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # const ax4a::Axis
    # ax4a_height::Observable{Float64}
    # ax4a_font::Observable{Any}
    # ax4a_fontsize::Observable{Float64}
    # ax4a_xlims::Observable{NTuple{2, Float64}}
    # ax4a_ylims::Observable{NTuple{2, Float64}}
    # ax4a_ion_positions::Observable{Vector{NTuple{2, Float64}}}
    # ax4a_ion_strings::Observable{Vector{String}}
    # ax4a_ion_colors::Observable{Vector{Symbol}}
    # ax4a_ion_offsets::Observable{Vector{NTuple{2, Float64}}}

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # AXIS 4a
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # const ax4b::Axis
    # ax4b_xlims::Observable{NTuple{2, Float64}}
    # ax4b_ylims::Observable{NTuple{2, Float64}}
    # ax4b_xticks_font::Observable{Any}
    # ax4b_xticks_fontsize::Observable{Float64}
    # ax4b_pixelwidth::Observable{Int}
    # ax4b_pixelheight::Observable{Int}
    # ax4b_yaxes_extra_stretchfactor::Observable{Float64}
    # ax4b_stretchfactor_changespeed_slow::Float64
    # ax4b_stretchfactor_changespeed_fast::Float64
    # ax4b_ms_ions::Observable{Vector{Float64}}
    # ax4b_ms_lineheights::Observable{Vector{Float64}}
    # ax4b_ms_linewidth::Observable{Float64}
    # ax4b_ms_linecolor::Observable{Symbol}
    # ax4b_ms_label_positions::Observable{Vector{NTuple{2, Float32}}}
    # ax4b_ms_label_texts::Observable{Vector{String}}
    # ax4b_ms_label_font::Observable{Any}
    # ax4b_ms_label_fontsize::Observable{Float64}
    # ax4b_ms_label_offset::Observable{NTuple{2, Float64}}
    # ax4b_slider_position::Observable{Vector{Float64}}
    # ax4b_slider_linewidth::Observable{Float64}
    # ax4b_slider_linecolor::Observable{Symbol}
    # ax4b_selected_ion_index::Observable{Union{Nothing, Int}}
    # ax4b_selected_ions::Observable{Set{SelectedIon}}
    # ax4b_selected_ion_positions::Observable{Vector{Float64}}
    # ax4b_selected_ion_colors::Observable{Vector{Symbol}}

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Fields related to the storage and selection of runsets.
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # runsets::Observable{Vector{RunSet}}
    # runset_index::Observable{Union{Int, Nothing}}
    # current_runset_cleanup::Observable{Bool}

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Fields related to the import of runsets
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # filereader::Union{Function, Nothing}

    ######################################################################################
    # Inner constructor method
    ######################################################################################
    function ExplorerData(
        ; reader::Union{Function, Nothing}=nothing,
        figure_size::NTuple{2, Int}=(1528, 750), 
        focus_on_show::Bool=true
        )

        ######################################################################################
        # Default values for various fields. These are later passed as values of keyworded
        # arguments to the constructor method.
        ######################################################################################

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Default values related to the window and its displayed figure in general.
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        figure_title_default = "JuChrom Explorer"
        figure_backgroundcolor_default = RGBA{Float32}(1.0f0,1.0f0,1.0f0,1.0f0)
        
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Default values related to Axes 1–4 equally
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ticklabelspace_left_default = 60
        # ticklabelspace_right_default = 60

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Default values related to Axis 1.
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ax1_font_default = "TeX Gyre Heros Makie Bold"
        # ax1_fontsize_default = 15
        # ax1_textcolor = :black
        # ax1_left_linecolor_default = :turquoise4
        # ax1_right_linecolor_default = :orange
        
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Default values related to Axis 2a.
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ax2a_font_default = "TeX Gyre Heros Makie Regular"
        # ax2a_fontsize_default = 15
        # ax2a_stored_data_marker_default = :vline
        # ax2a_stored_data_markersize_default = 20
        # ax2b_stored_data_color_default = :red

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Default values not related or related equally to the left and right y-2b axes.
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ax2b_xaxis_datatype_default = :time
        # ax2b_xaxis_timeunit_default = u"minute"
        # ax2b_xticks_font_default = "TeX Gyre Heros Makie Regular"
        # ax2b_xticks_fontsize_default = 14
        # ax2b_yticks_font_default = "TeX Gyre Heros Makie Regular"
        # ax2b_yticks_fontsize_default = 14
        # ax2b_yaxes_extra_stretchfactor_default = 1.05
        # ax2b_stretchfactor_changespeed_slow_default = 1.05
        # ax2b_stretchfactor_changespeed_fast_default = 1.5
        # ax2b_linewidth_for_current_runset_default = 1.5
        # ax2b_linewidth_for_non_current_runset_default = 0.5
        # ax2b_crosshair_color_default = :grey
        # ax2b_crosshair_linewidth_default = 0.5
        # ax2b_crosshair_label_font_default = "TeX Gyre Heros Makie Regular"
        # ax2b_crosshair_label_fontsize_default = 12
        # ax2b_gcmsdata_scantimes_add_linewidth_default = 1.5
        # ax2b_gcmsdata_scantimes_add_color_default = :grey5
        # ax2b_gcmsdata_scantimes_subtract_linewidth_default = 1.5
        # ax2b_gcmsdata_scantimes_subtract_color_default = :orangered3

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Default values related to axis 3a.
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ax3a_font_default = "TeX Gyre Heros Makie Regular"
        # ax3a_fontsize_default = 14
        # ax3b_ms_linewidth_default = 1
        # ax3b_ms_linecolor_default = :orangered3
        # ax3b_ms_label_font_default = "TeX Gyre Heros Makie Regular"
        # ax3b_ms_label_fontsize_default = 12
        # ax3b_ms_label_offset_default = (0.0, 5.0)

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Default values related to axis 3b.
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ax3b_ion_chromatogram_linewidth_default = 1
        # ax3b_ion_chromatogram_color_default = :grey
        # ax3b_crosshair_color_default = :grey
        # ax3b_crosshair_linewidth_default = 0.5

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Default values related to axis 4a.
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ax4a_font_default = "TeX Gyre Heros Makie Regular"
        # ax4a_fontsize_default = 14

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Default values related to axis 4b.
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ax4b_xticks_font_default = "TeX Gyre Heros Makie Regular"
        # ax4b_xticks_fontsize_default = 15
        # ax4b_stretchfactor_changespeed_slow_default = 1.05
        # ax4b_stretchfactor_changespeed_fast_default = 1.5
        # ax4b_ms_linewidth_default = 1
        # ax4b_ms_linecolor_default = :orangered3
        # ax4b_ms_label_font_default = "TeX Gyre Heros Makie Regular"
        # ax4b_ms_label_fontsize_default = 12
        # ax4b_ms_label_offset_default = (0.0, 5.0)
        # ax4b_slider_linewidth_default = 0.5
        # ax4b_slider_linecolor_default = :grey

        ######################################################################################
        # Perform all necessary preparations to create a new ExplorerData instance.
        ######################################################################################

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Create a interative GLMakie figure
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        first(figure_size) > 0 || throw(
            ArgumentError("figure width must be a positive integer larger than zero"))
        last(figure_size) > 0 || throw(
            ArgumentError("figure height must be a positive integer larger than zero"))

        GLMakie.activate!(; title=figure_title_default, focus_on_show=focus_on_show)
        fig = Figure(; backgroundcolor=figure_backgroundcolor_default, size=figure_size)

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Create eight axes in four rows:
        # subfig 1: ax1
        # subfig 1: ax2a
        #           ax2b_left and ax2b_right, sharing the same x-axis
        # subfig 3: ax3a
        #           ax3b
        # subfig 4: ax4a
        #           ax4b
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # ticklabelspace_left = Observable{Float64}(ticklabelspace_left_default)
        # ticklabelspace_right = Observable{Float64}(ticklabelspace_right_default)

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # SUBFIG 1: Axis 1
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # subfig1 = fig[1, 1] = GridLayout()
        
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # AXIS 1
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # The height of the axis is derived from the size of the text it displays. Therefore,
        # the font and its size must be handled first.
        # ax1_font = Observable(ax1_font_default)
        # ax1_fontsize = Observable{Float64}(ax1_fontsize_default)

        # Get the height of the axis based on the size of the text it displays.
        # ax1_height = fontheight(ax1_font, ax1_fontsize)

        # Create axis with predefined height.
        # ax1 = Axis(subfig1[1, 1], height=ax1_height)

        # Explicitly set axis limits to default values.
        # limits!(ax1, 0, 10, 0, 10)

        # Disable all mouse interactions (e.g., zooming, scrolling).
        # deregister_default_mouse_interactions!(ax1)

        # Axis does not need decorations and splines.
        # hidedecorations!(ax1)
        # hidespines!(ax1)

        # Create observales and a text plot to store and display the name of the runset.
        # ax1_run_name = Observable{String}("")
        # ax1_textcolor = Observable{Symbol}(ax1_textcolor)
        # text!(ax1,
        #     Point2f([0.5, 0.5]),
        #     text=ax1_run_name,            # ::Observable
        #     font=ax1_font,                # ::Observable
        #     fontsize=ax1_fontsize,        # ::Observable
        #     color=ax1_textcolor,          # ::Observable
        #     align=(:center, :center),
        #     space=:relative)

        # Create observales and two poly plots to store and display the default colors of the lines
        # displayed on the left and right y-2b axes.
        # ax1_left_linecolor = Observable{Symbol}(ax1_left_linecolor_default)
        # ax1_right_linecolor = Observable{Symbol}(ax1_right_linecolor_default)
        # poly!(ax1, leftcolorbox(ax1_height), color=ax1_left_linecolor, space = :pixel)
        # poly!(ax1, rightcolorbox(ax1, ax1_height), color=ax1_right_linecolor, space = :pixel)

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # SUBFIG 2: Axis 2a, Axis2b_left, Axis2b_right
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # subfig2 = fig[2, 1] = GridLayout()

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # AXIS 2a
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # The height of the axis is derived from the size of the text it displays. Therefore,
        # the font and its size must be handled first.
        # ax2a_font = Observable(ax2a_font_default)
        # ax2a_fontsize = Observable{Float64}(ax2a_fontsize_default)
        
        # Get the height of the axis based on the size of the text it displays.
        # ax2a_height = fontheight(ax2a_font, ax2a_fontsize)
        
        # Create axis with predefined height.
        # ax2a = Axis(subfig2[1, 1], height=ax2a_height)

        # Explicitly set axis limits to default values.
        # limits!(ax2a, 0, 10, 0, 10)

        # Disable all mouse interactions (e.g., zooming, scrolling).
        # deregister_default_mouse_interactions!(ax2a)

        # Axis does not need decorations and splines.
        # hidedecorations!(ax2a)
        # hidespines!(ax2a)
        
        # ax2a_stored_data_pos_x_y = Observable{Vector{NTuple{2, Float64}}}([])
        # ax2a_stored_data_marker = Observable(ax2a_stored_data_marker_default)
        # ax2a_stored_data_markersize = Observable(ax2a_stored_data_markersize_default)
        # ax2a_stored_data_color = Observable(ax2b_stored_data_color_default)

        # ax2a_stored_data_plot = scatter!(ax2a,
        #     ax2a_stored_data_pos_x_y,
        #     marker=ax2a_stored_data_marker,
        #     markersize=ax2a_stored_data_markersize,
        #     color=ax2a_stored_data_color)     

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # AXIS 2a left
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ax2b_left = Axis(subfig2[2, 1])

        # Disable all mouse interactions (e.g., zooming, scrolling).
        # deregister_default_mouse_interactions!(ax2b_left)

        # Hide the axis decorations, since we do not know if the axis will contain data to
        # display.
        # hidedecorations!(ax2b_left)

        # Label the y-axis as it will always display abundance data.
        # ax2b_left.ylabel = "Abundance"
        # on(ticklabelspace_left, update=true, priority=1) do ticklabelspace
        #     ax2b_left.yaxis.attributes.ticklabelspace[] = ticklabelspace
        # end
        
        # Create observables to store whether the y-axis is needed, what its current limits are,
        # and by what and by what factor the axis should be stretched.
        # ax2b_left_required = Observable{Bool}(false)
        # ax2b_left_ylims = Observable{NTuple{2, Float64}}((0, 10))
        # on(ax2b_left_ylims, update=true, priority=1) do ylims_left
        #     ylims!(ax2b_left,  ylims_left...)
        # end
        # ax2b_left_yaxis_stretchfactor = Observable{Float64}(1.0)

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # AXIS 2a right
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ax2b_right = Axis(subfig2[2, 1], yaxisposition=:right)

        # Disable all mouse interactions (e.g., zooming, scrolling).
        # deregister_default_mouse_interactions!(ax2b_right)

        # Hide the axis decorations, since we do not know if the axis will contain data to
        # display.
        # hidedecorations!(ax2b_right)

        # Hide the splines of the right axis, since all necessary splines are already displayed
        # by the left axis.
        # hidespines!(ax2b_right)
    
        # Label the y-axis as it will always display abundance data.
        # ax2b_right.ylabel = "Abundance"

        # on(ticklabelspace_right, update=true, priority=1) do ticklabelspace
        #     ax2b_right.yaxis.attributes.ticklabelspace[] = ticklabelspace
        # end

        # Create observables to store whether the y-axis is needed, what its current limits are,
        # and by what and by what factor the axis should be stretched.
        # ax2b_right_required = Observable{Bool}(false)
        # ax2b_right_ylims = Observable{NTuple{2, Float64}}((0, 10))
        # on(ax2b_right_ylims, update=true, priority=1) do ylims_right
        #     ylims!(ax2b_right,  ylims_right...)
        # end
        # ax2b_right_yaxis_stretchfactor = Observable{Float64}(1.0)

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Ax2b – data or actions that are not specific to either the left or right axis
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # 2b left axis and 2b right axis share a common x axis
        # linkxaxes!(ax2b_left, ax2b_right)

        # Create observables to store whether the x-axis is needed and to characterize
        # the ax2 in general.
        # ax2b_xaxis_required = Observable{Bool}(false)
        # ax2b_xaxis_datatype = Observable{Symbol}(ax2b_xaxis_datatype_default)
        # ax2b_xaxis_minintervalsize = Observable{Float64}(10^-5)
        # ax2b_xaxis_timeunit = Observable{Unitful.TimeUnits}(ax2b_xaxis_timeunit_default)
        # ax2b_xlims = Observable{NTuple{2, Float64}}((0, 10))
        # on(ax2b_xlims, update=true, priority=1) do xlims
        #     xlims!(ax2b_left,  xlims...)
        #     xlims!(ax2b_right,  xlims...)
        # end
        
        # ax2a has the same x limits as ax2b
        # on(ax2b_xlims, update=true, priority=1) do xlims
        #     xlims!(ax2a,  xlims...)
        # end

        # ax2b_xticks_font = Observable{Any}(ax2b_xticks_font_default)
        # ax2b_xticks_fontsize = Observable{Float64}(ax2b_xticks_fontsize_default)
        # onany(ax2b_xticks_font, ax2b_xticks_fontsize, update=true, priority=1) do font, fontsize
        #     ax2b_left.xticklabelfont = font
        #     ax2b_left.xticklabelsize = fontsize
        # end

        # ax2b_yticks_font = Observable{Any}(ax2b_yticks_font_default)
        # ax2b_yticks_fontsize = Observable{Float64}(ax2b_yticks_fontsize_default)
        # onany(ax2b_yticks_font, ax2b_yticks_fontsize, update=true, priority=1) do font, fontsize
        #     ax2b_left.yticklabelfont = font
        #     ax2b_left.yticklabelsize = fontsize
        #     ax2b_right.yticklabelfont = font
        #     ax2b_right.yticklabelsize = fontsize
        # end

        # ax2b_pixelwidth = size_in_px_x(ax2b_left)
        # ax2b_pixelheight = size_in_px_y(ax2b_left)

        # ax2b_yaxes_extra_stretchfactor = Observable{Float64}(ax2b_yaxes_extra_stretchfactor_default)
        # ax2b_stretchfactor_changespeed_slow = ax2b_stretchfactor_changespeed_slow_default
        # ax2b_stretchfactor_changespeed_fast = ax2b_stretchfactor_changespeed_fast_default

        # ax2b_linewidth_for_current_runset = Observable{Float64}(ax2b_linewidth_for_current_runset_default)
        # ax2b_linewidth_for_non_current_runset = Observable{Float64}(ax2b_linewidth_for_non_current_runset_default)

        # Create observables and a vertical and a horizonal line plot to store data on and display
        # a crosshair
        # ax2b_crosshair_pos_x = Observable{Vector{Float64}}([])
        # ax2b_crosshair_pos_y = Observable{Vector{Float64}}([])
        # ax2b_crosshair_color = Observable{Symbol}(ax2b_crosshair_color_default)
        # ax2b_crosshair_linewidth = Observable{Float64}(ax2b_crosshair_linewidth_default)

        # ax2b_crosshair_plot_x = vlines!(ax2b_left,
        #     ax2b_crosshair_pos_x,                   # ::Observable
        #     color=ax2b_crosshair_color,             # ::Observable
        #     linewidth=ax2b_crosshair_linewidth)     # ::Observable
        
        # ax2b_crosshair_plot_y = hlines!(ax2b_left,
        #     ax2b_crosshair_pos_y,                   # ::Observable
        #     color=ax2b_crosshair_color,             # ::Observable
        #     linewidth=ax2b_crosshair_linewidth)     # ::Observable
        
        # translate!(ax2b_crosshair_plot_x, 0, 0, 10^4)  # Keep crosshair in the foreground!
        # translate!(ax2b_crosshair_plot_y, 0, 0, 10^4)  # Keep crosshair in the foreground!

        # ax2b_dragstart = Vector{Float64}([])
        # ax2b_left_prevlimits = Vector{NTuple{4, Float64}}()
        # ax2b_right_prevlimits = Vector{NTuple{4, Float64}}()

        # Create observables and a text plot to display information on the crosshair
        # ax2b_crosshair_label_font = Observable(ax2b_crosshair_label_font_default)
        # ax2b_crosshair_label_fontsize = Observable{Float64}(ax2b_crosshair_label_fontsize_default)
        # ax2b_crosshair_label_positions = Observable{Vector{NTuple{2, Float64}}}([])
        # ax2b_crosshair_label_texts = Observable{Vector{String}}([])
        # ax2b_crosshair_label_alignments = Observable{Vector{NTuple{2, Symbol}}}([])
        # ax2b_crosshair_label_offsets = Observable{Vector{NTuple{2, Float64}}}([])
        # ax2b_crosshair_label_visible = Observable{Bool}(false)

        # ax2b_crosshair_label_plot = text!(ax2b_left,
        #     ax2b_crosshair_label_positions,          # ::Observable
        #     text=ax2b_crosshair_label_texts,         # ::Observable
        #     align=ax2b_crosshair_label_alignments,   # ::Observable
        #     offset=ax2b_crosshair_label_offsets,     # ::Observable
        #     font=ax2b_crosshair_label_font,          # ::Observable
        #     fontsize=ax2b_crosshair_label_fontsize,  # ::Observable
        #     color=:black,
        #     glowcolor=:white,
        #     glowwidth=3)
        
        # translate!(ax2b_crosshair_label_plot, 0, 0, 10^4)  # Keep crosshair labels in the foreground!

        # Create observables and two vertical line plots to store data about which scans' intensities
        # to add and which scans' intensities to subtract when extracting a mass spectrum.
        # ax2b_gcmsdata_scanindexselect_start = Observable{Union{Nothing, Int}}(nothing)
        # ax2b_gcmsdata_scanindexselect_stop = Observable{Union{Nothing, Nothing}}(nothing)
        # ax2b_gcmsdata_scanindexselect_current_range = Observable{Union{Nothing, UnitRange{Int64}}}(nothing)
        # ax2b_gcmsdata_scanindexselect_mode = Observable{Union{Nothing, Symbol}}(nothing)

        # ax2b_gcmsdata_scanindexselect_add = Observable{Set{Int}}(Set([]))
        # ax2b_gcmsdata_scanindexselect_subtract = Observable{Set{Int}}(Set([]))

        # ax2b_gcmsdata_scantimes_add = Observable{Vector{Float64}}([])
        # ax2b_gcmsdata_scantimes_add_linewidth = Observable{Float64}(ax2b_gcmsdata_scantimes_add_linewidth_default)
        # ax2b_gcmsdata_scantimes_add_color = Observable{Symbol}(ax2b_gcmsdata_scantimes_add_color_default)

        # vlines!(ax2b_left,
        #     ax2b_gcmsdata_scantimes_add,                      # ::Observable
        #     linewidth=ax2b_gcmsdata_scantimes_add_linewidth,  # ::Observable
        #     color=ax2b_gcmsdata_scantimes_add_color)          # ::Observable

        # ax2b_gcmsdata_scantimes_subtract = Observable{Vector{Float64}}([])
        # ax2b_gcmsdata_scantimes_subtract_linewidth = Observable{Float64}(ax2b_gcmsdata_scantimes_subtract_linewidth_default)
        # ax2b_gcmsdata_scantimes_subtract_color = Observable{Symbol}(ax2b_gcmsdata_scantimes_subtract_color_default)

        # vlines!(ax2b_left,
        #     ax2b_gcmsdata_scantimes_subtract,                      # ::Observable
        #     linewidth=ax2b_gcmsdata_scantimes_subtract_linewidth,  # ::Observable
        #     color=ax2b_gcmsdata_scantimes_subtract_color)          # ::Observable

        # No gap between ax2a and ax2b_left/_right1
        # rowgap!(subfig2, 0)

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # SUBFIG 3: Axis 3a, Axis3b
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # subfig3 = fig[3, 1] = GridLayout()

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # AXIS 3a
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # The height of the axis is derived from the size of the text it displays. Therefore,
        # the font and its size must be handled first.
        # ax3a_font = Observable(ax3a_font_default)
        # ax3a_fontsize = Observable{Float64}(ax3a_fontsize_default)
    
        # Get the height of the axis based on the size of the text it displays.
        # ax3a_height = fontheight(ax3a_font, ax3a_fontsize)
        
        # Create axis with predefined height.
        # ax3a = Axis(subfig3[1, 1], height=ax3a_height)
        
        # Disable all mouse interactions (e.g., zooming, scrolling).
        # deregister_default_mouse_interactions!(ax3a)
        
        # Axis does not need decorations and splines.
        # hidedecorations!(ax3a)
        # hidespines!(ax3a)
        
        # Create observables to characterize ax3a in general.
        # ax3a_xlims = Observable{NTuple{2, Float64}}((0, 100))
        # on(ax3a_xlims, update=true, priority=1) do xlims
        #     xlims!(ax3a, xlims...)
        # end
        # ax3a_ylims = Observable{NTuple{2, Float64}}((0, 10))
        # on(ax3a_ylims, update=true, priority=1) do ylims
        #     ylims!(ax3a, ylims...)
        # end

        # Create observables and a text plot to store and name the ions whose chromatograms are displayed.
        # ax3a_extracted_ion_strings = Observable{Vector{Makie.RichText}}([])
        # ax3a_extracted_ion_string_positions = Observable{Vector{NTuple{2, Float64}}}([])

        # ax3a_extracted_ion_string_plot = text!(ax3a,
        #     ax3a_extracted_ion_string_positions,  # ::Observable
        #     text=ax3a_extracted_ion_strings,      # ::Observable
        #     font=ax3a_font,                       # ::Observable
        #     fontsize=ax3a_fontsize,               # ::Observable
        #     offset=(0.0, 1.0),
        #     align=(:left, :bottom))

        # Create observables and a text plot to store and name the ions whose chromatograms are displayed.
        # ax3a_ion_intensity_extrema_strings = Observable{Vector{String}}([])
        # ax3a_ion_intensity_extrema_string_positions = Observable{Vector{NTuple{2, Float64}}}([])

        # text!(ax3a,
        #     ax3a_ion_intensity_extrema_string_positions,  # ::Observable
        #     text=ax3a_ion_intensity_extrema_strings,      # ::Observable
        #     font=ax3a_font,                               # ::Observable
        #     fontsize=ax3a_fontsize,                       # ::Observable
        #     offset=(0.0, 1.0),
        #     align=(:right, :bottom))

        


        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # AXIS 3b
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ax3b = Axis(subfig3[2, 1])

        # Disable all mouse interactions (e.g., zooming, scrolling).
        # deregister_default_mouse_interactions!(ax3b)
        
        # Hide the axis decorations since we do not know what kind of data will be displayed.
        # hidedecorations!(ax3b)

        # Set a default tick label space for the y-axis.
        # on(ticklabelspace_left, update=true, priority=1) do ticklabelspace
        #     ax3b.yaxis.attributes.ticklabelspace[] = ticklabelspace
        # end

        # Create an observables to monitor what type of data axis 3b should display.
        # ax3b_mode = Observable{Symbol}(:undefined)
        
        # Create observables to characterize ax3b in general.
        # ax3b_xlims = Observable{NTuple{2, Float64}}((0, 10))
        # on(ax3b_xlims, update=true, priority=1) do xlims
        #     xlims!(ax3b, xlims...)
        # end
        # ax3b_ylims = Observable{NTuple{2, Float64}}((0, 10))
        # on(ax3b_ylims, update=true, priority=1) do ylims
        #     ylims!(ax3b, ylims...)
        # end
        # ax3b_pixelwidth = size_in_px_x(ax3b)
        # ax3b_pixelheight = size_in_px_y(ax3b)

        # Create observables, a vline plot, and a text plot to store data for and enable a
        # mass spectrum live view.
        # ax3b_ms_ions = Observable{Vector{Float64}}([])
        # ax3b_ms_lineheights = Observable{Vector{Float64}}([])
        # # ax3b_ms_linewidth = Observable{Float64}(ax3b_ms_linewidth_default)
        # ax3b_ms_linecolor = Observable{Symbol}(ax3b_ms_linecolor_default) 
        # ax3b_ms_label_positions = Observable{Vector{NTuple{2, Float32}}}([])
        # ax3b_ms_label_texts = Observable{Vector{String}}([])
        # ax3b_ms_label_font = Observable{Any}(ax3b_ms_label_font_default)
        # ax3b_ms_label_fontsize = Observable{Float64}(ax3b_ms_label_fontsize_default)
        # ax3b_ms_label_offset = Observable{NTuple{2, Int}}(ax3b_ms_label_offset_default) 

        # ax3b MS is always formatted as the MS of ax4b
        # ax4b_ms_linewidth = Observable{Float64}(ax4b_ms_linewidth_default)
        # ax4b_ms_label_font = Observable{Any}(ax4b_ms_label_font_default)
        # ax4b_ms_label_fontsize = Observable{Float64}(ax4b_ms_label_fontsize_default)
        # ax4b_ms_label_offset = Observable{NTuple{2, Int}}(ax4b_ms_label_offset_default) 


        # vlines!(ax3b,
        #     ax3b_ms_ions,                     # ::Observable
        #     ymax=ax3b_ms_lineheights,         # ::Observable
        #     linewidth=ax4b_ms_linewidth,      # ::Observable
        #     color=ax3b_ms_linecolor)          # ::Observable

        # text!(ax3b,
        #     ax3b_ms_label_positions,          # ::Observable
        #     text=ax3b_ms_label_texts,         # ::Observable
        #     font=ax4b_ms_label_font,          # ::Observable
        #     fontsize=ax4b_ms_label_fontsize,  # ::Observable
        #     offset=ax4b_ms_label_offset,      # ::Observable
        #     align=(:center, :bottom))

        # Create observables and a lines plot to store data for and enable an
        # ion chromatogram live view.
        # ax3b_ion_chromatogram_scantimes = Observable{Vector{Float64}}([0])    # Note: "lines!" with an empty vector for x causes Makie to crash.
        # ax3b_ion_chromatogram_intensities = Observable{Vector{Float64}}([0])  # Note: "lines!" with an empty vector for x causes Makie to crash.
        # ax3b_ion_chromatogram_linewidth = Observable{Float64}(ax3b_ion_chromatogram_linewidth_default)
        # ax3b_ion_chromatogram_color = Observable{Symbol}(ax3b_ion_chromatogram_color_default)

        # ax3b_ion_chromatogram = lines!(ax3b,
        #     ax3b_ion_chromatogram_scantimes,
        #     ax3b_ion_chromatogram_intensities,
        #     linewidth=ax3b_ion_chromatogram_linewidth,
        #     color=ax3b_ion_chromatogram_color)
        
        # translate!(ax3b_ion_chromatogram, 0, 0, 10^4)  # Keep the live view ion chromtogram in the foreground!


        # Create observables for vertical crosshair
        # ax3b_crosshair_pos_x = Observable{Vector{Float64}}([])
        # ax3b_crosshair_color = Observable{Symbol}(ax3b_crosshair_color_default)
        # ax3b_crosshair_linewidth = Observable{Float64}(ax3b_crosshair_linewidth_default)

        # ax3b_crosshair_plot_x = vlines!(ax3b,
        #     ax3b_crosshair_pos_x,                   # ::Observable
        #     color=ax3b_crosshair_color,             # ::Observable
        #     linewidth=ax3b_crosshair_linewidth)     # ::Observable
        
        # translate!(ax3b_crosshair_plot_x, 0, 0, 10^4)  # Keep crosshair in the foreground!

        # No gap between ax3a and ax3b
        # rowgap!(subfig3, 0)

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # SUBFIG 4: Axis 4a, Axis4b
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # subfig4 = fig[4, 1] = GridLayout()
        
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # AXIS 4a
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # The height of the axis is derived from the size of the text it displays. Therefore,
        # the font and its size must be handled first.
        # ax4a_font = Observable(ax4a_font_default)
        # ax4a_fontsize = Observable{Float64}(ax4a_fontsize_default)
        
        # Get the height of the axis based on the size of the text it displays.
        # ax4a_height = fontheight(ax4a_font, ax4a_fontsize)
        
        # Create axis with predefined height.
        # ax4a = Axis(subfig4[1, 1], height=ax4a_height)
        
        # Disable all mouse interactions (e.g., zooming, scrolling).
        # deregister_default_mouse_interactions!(ax4a)
        
        # Axis does not need decorations and splines.
        # hidedecorations!(ax4a)
        # hidespines!(ax4a)

        # Create observables to characterize ax4a in general.
        # ax4a_xlims = Observable{NTuple{2, Float64}}((0, 10))
        # on(ax4a_xlims, update=true, priority=1) do xlims
        #     xlims!(ax4a, xlims...)
        # end
        # ax4a_ylims = Observable{NTuple{2, Float64}}((0, 10))
        # on(ax4a_ylims, update=true, priority=1) do ylims
        #     ylims!(ax4a, ylims...)
        # end

        # Create observables and a text plot to store data for and name the ion
        # on which the cursor is pointing in axis 4b.
        # ax4a_ion_positions = Observable{Vector{NTuple{2, Float64}}}([])
        # ax4a_ion_strings = Observable{Vector{String}}([])
        # ax4a_ion_colors = Observable{Vector{Symbol}}([])
        # ax4a_ion_offsets = Observable{Vector{NTuple{2, Float64}}}([])

        # text!(ax4a,
        #     ax4a_ion_positions,       # ::Observable
        #     text=ax4a_ion_strings,    # ::Observable
        #     font=ax4a_font,           # ::Observable
        #     fontsize=ax4a_fontsize,   # ::Observable
        #     color=ax4a_ion_colors,    # ::Observable
        #     offset=ax4a_ion_offsets,  # ::Observable
        #     align=(:center, :bottom))  

        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # AXIS 4b
        # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # ax4b = Axis(subfig4[2, 1])
        
        # Disable all mouse interactions (e.g., zooming, scrolling).
        # deregister_default_mouse_interactions!(ax4b)
        
        # Hide the axis decorations, since we do not know if the axis will contain data to
        # display.
        # hidedecorations!(ax4b)
        
        # Label the x-axis and y-axis because they will always display the same type of data.
        # ax4b.xlabel = rich(rich("m", font = :italic), rich("/", font = :regular), rich("z", font = :italic))
        # ax4b.ylabel = "Abundance [%]"
        # on(ticklabelspace_left, update=true, priority=1) do ticklabelspace
        #     ax4b.yaxis.attributes.ticklabelspace[] = ticklabelspace
        # end

        # Create observables to characterize ax4b in general.
        # ax4b_xlims = Observable{NTuple{2, Float64}}((0, 10))
        # on(ax4b_xlims, update=true, priority=1) do xlims
        #     xlims!(ax4b, xlims...)
        # end
        # ax4b_ylims = Observable{NTuple{2, Float64}}((0, 10))
        # on(ax4b_ylims, update=true, priority=1) do ylims
        #     ylims!(ax4b, ylims...)
        # end

        # ax4b_xticks_font = Observable{Any}(ax4b_xticks_font_default)
        # ax4b_xticks_fontsize = Observable{Float64}(ax4b_xticks_fontsize_default)
        # ax4b_pixelwidth = size_in_px_x(ax4b)
        # ax4b_pixelheight = size_in_px_y(ax4b)

        # Observables that store the factor by which the y-axis is stretched by the user,
        # and the speed at which the stretch factor can be changed using keyboard shortcuts.
        # ax4b_yaxes_extra_stretchfactor = Observable{Float64}(1.0)
        # ax4b_stretchfactor_changespeed_slow = ax4b_stretchfactor_changespeed_slow_default
        # ax4b_stretchfactor_changespeed_fast = ax4b_stretchfactor_changespeed_fast_default


        # Create observables, a vline plot, and a text plot to store data for and display a
        # mass spectrum extracted from scans selected in ax2b.
        # ax4b_ms_ions = Observable{Vector{Float64}}([])
        # ax4b_ms_lineheights = Observable{Vector{Float64}}([])
        # ax4b_ms_linewidth = Observable{Float64}(ax4b_ms_linewidth_default)
        # ax4b_ms_linecolor = Observable{Symbol}(ax4b_ms_linecolor_default) 
        # ax4b_ms_label_positions = Observable{Vector{NTuple{2, Float32}}}([])
        # ax4b_ms_label_texts = Observable{Vector{String}}([])

        # vlines!(ax4b,
        #     ax4b_ms_ions,                     # ::Observable
        #     ymax=ax4b_ms_lineheights,         # ::Observable
        #     linewidth=ax4b_ms_linewidth,      # ::Observable
        #     color=ax4b_ms_linecolor)          # ::Observable
        # text!(ax4b, 
        #     ax4b_ms_label_positions,          # ::Observable
        #     text=ax4b_ms_label_texts,         # ::Observable
        #     font=ax4b_ms_label_font,          # ::Observable
        #     fontsize=ax4b_ms_label_fontsize,  # ::Observable
        #     offset=ax4b_ms_label_offset,      # ::Observable
        #     align=(:center, :bottom))

        # Create observables and a vlines plot to store data for and display a slider
        # to select an ion whose chromtogram is displayed in live view on axis 3b.
        # ax4b_slider_position = Observable{Vector{Float64}}([])
        # ax4b_slider_linewidth = Observable{Float64}(ax4b_slider_linewidth_default)
        # ax4b_slider_linecolor = Observable{Symbol}(ax4b_slider_linecolor_default)

        # ax4b_slider_plot = vlines!(ax4b,
        #     ax4b_slider_position,             # ::Observable
        #     linewidth=ax4b_slider_linewidth,  # ::Observable
        #     color=ax4b_slider_linecolor)      # ::Observable
        # translate!(ax4b_slider_plot, 0, 0, 10^4)

        # ax4b_selected_ion_index = Observable{Union{Nothing, Int}}(nothing)
        
        # Create observables to store information about which chromtograms of selected
        # ions to display on axis 3b.
        # ax4b_selected_ions = Observable{Set{SelectedIon}}(Set([]))
        # ax4b_selected_ion_positions = Observable{Vector{Float64}}([])
        # ax4b_selected_ion_colors = Observable{Vector{Symbol}}([])

        # No gap between ax4a and ax4b
        # rowgap!(subfig4, 0)

        # Create observables to store and manage runsets.
        # runsets = Observable{Vector{RunSet}}([])
        # runset_index = Observable{Union{Int, Nothing}}(nothing)
        # current_runset_cleanup = Observable{Bool}(false)

        # Store any provided default filereader.
        # filereader = reader

        # Display the interactive figure in a window.
        window = GLMakie.to_native(display(fig))

        ######################################################################################
        # Create a new Explorer data instance
        ######################################################################################
        new(
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # Fields related to the window and its displayed figure in general.
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            fig,
            Observable(figure_backgroundcolor_default),
            window # ,

            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # AXIS 1–4
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ticklabelspace_left,
            # ticklabelspace_right,

            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # AXIS 1
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ax1,
            # ax1_height,
            # ax1_font,
            # ax1_fontsize,
            # ax1_textcolor,
            # ax1_run_name,
            # ax1_left_linecolor,
            # ax1_right_linecolor,

            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # AXIS 2a
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ax2a,
            # ax2a_height,
            # ax2a_font,
            # ax2a_fontsize,
            # ax2a_stored_data_pos_x_y,
            # ax2a_stored_data_marker,
            # ax2a_stored_data_markersize,
            # ax2a_stored_data_color,
            # ax2a_stored_data_plot,

            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # AXIS 2b - general
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
            # ax2b_xaxis_required,
            # ax2b_xaxis_datatype,
            # ax2b_xaxis_minintervalsize,
            # ax2b_xaxis_timeunit,
            # ax2b_xlims,
            # ax2b_xticks_font,
            # ax2b_xticks_fontsize,
            # ax2b_yticks_font,
            # ax2b_yticks_fontsize,
            # ax2b_pixelwidth,
            # ax2b_pixelheight,
            # ax2b_yaxes_extra_stretchfactor,
            # ax2b_stretchfactor_changespeed_slow,
            # ax2b_stretchfactor_changespeed_fast,
            # ax2b_linewidth_for_current_runset,
            # ax2b_linewidth_for_non_current_runset,
            # ax2b_crosshair_pos_x,
            # ax2b_crosshair_pos_y,
            # ax2b_crosshair_color,
            # ax2b_crosshair_linewidth,
            # ax2b_dragstart,
            # ax2b_crosshair_label_font,
            # ax2b_crosshair_label_fontsize,
            # ax2b_crosshair_label_positions,
            # ax2b_crosshair_label_texts,
            # ax2b_crosshair_label_alignments,
            # ax2b_crosshair_label_offsets,
            # ax2b_crosshair_label_visible,
            # ax2b_gcmsdata_scanindexselect_start,
            # ax2b_gcmsdata_scanindexselect_stop,
            # ax2b_gcmsdata_scanindexselect_current_range,
            # ax2b_gcmsdata_scanindexselect_mode,
            # ax2b_gcmsdata_scanindexselect_add,
            # ax2b_gcmsdata_scanindexselect_subtract,
            # ax2b_gcmsdata_scantimes_add,
            # ax2b_gcmsdata_scantimes_add_linewidth,
            # ax2b_gcmsdata_scantimes_add_color,
            # ax2b_gcmsdata_scantimes_subtract,
            # ax2b_gcmsdata_scantimes_subtract_linewidth,
            # ax2b_gcmsdata_scantimes_subtract_color,

            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # AXIS 2b - left y axis
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ax2b_left,
            # ax2b_left_required,
            # ax2b_left_ylims,
            # ax2b_left_yaxis_stretchfactor,
            # ax2b_left_prevlimits,

            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # AXIS 2b - right y axis
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ax2b_right,
            # ax2b_right_required,
            # ax2b_right_ylims,
            # ax2b_right_yaxis_stretchfactor,
            # ax2b_right_prevlimits,

            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # AXIS 3a
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ax3a,
            # ax3a_height,
            # ax3a_font,
            # ax3a_fontsize,
            # ax3a_xlims,
            # ax3a_ylims,        
            # ax3a_extracted_ion_strings,
            # ax3a_extracted_ion_string_positions,
            # ax3a_extracted_ion_string_plot,
            # ax3a_ion_intensity_extrema_strings,
            # ax3a_ion_intensity_extrema_string_positions,
        
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # AXIS 3b
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ax3b,
            # ax3b_xlims,
            # ax3b_ylims,
            # ax3b_pixelwidth,
            # ax3b_pixelheight,
            # ax3b_mode,
            # ax3b_ms_ions,
            # ax3b_ms_lineheights,
            # # ax3b_ms_linewidth,
            # ax3b_ms_linecolor,
            # ax3b_ms_label_positions,
            # ax3b_ms_label_texts,
            # # ax3b_ms_label_font,
            # # ax3b_ms_label_fontsize,
            # # ax3b_ms_label_offset,
            # ax3b_ion_chromatogram_scantimes,
            # ax3b_ion_chromatogram_intensities,
            # ax3b_ion_chromatogram_linewidth,
            # ax3b_ion_chromatogram_color,
            # ax3b_crosshair_pos_x,
            # ax3b_crosshair_color,
            # ax3b_crosshair_linewidth,

            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # AXIS 4a
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ax4a,
            # ax4a_height,
            # ax4a_font,
            # ax4a_fontsize,
            # ax4a_xlims,
            # ax4a_ylims,    
            # ax4a_ion_positions,
            # ax4a_ion_strings,
            # ax4a_ion_colors,
            # ax4a_ion_offsets,

            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # AXIS 4b
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ax4b,
            # ax4b_xlims,
            # ax4b_ylims,
            # ax4b_xticks_font,
            # ax4b_xticks_fontsize,
            # ax4b_pixelwidth,
            # ax4b_pixelheight,
            # ax4b_yaxes_extra_stretchfactor,
            # ax4b_stretchfactor_changespeed_slow,
            # ax4b_stretchfactor_changespeed_fast,
            # ax4b_ms_ions,
            # ax4b_ms_lineheights,
            # ax4b_ms_linewidth,
            # ax4b_ms_linecolor,
            # ax4b_ms_label_positions,
            # ax4b_ms_label_texts,
            # ax4b_ms_label_font,
            # ax4b_ms_label_fontsize,
            # ax4b_ms_label_offset,
            # ax4b_slider_position,
            # ax4b_slider_linewidth,
            # ax4b_slider_linecolor,
            # ax4b_selected_ion_index,
            # ax4b_selected_ions,
            # ax4b_selected_ion_positions,
            # ax4b_selected_ion_colors,

            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # Fields related to the storage and selection of runsets.
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # runsets,
            # runset_index,
            # current_runset_cleanup,

            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # Fields related to the import of runsets
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            # filereader
        )
    end  # of inner constructor method
end  # of mutable struct

######################################################################################
# Getters and setters of ExplorerData fields.
######################################################################################

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Fields related to the window and its displayed figure in general.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fig(e::ExplorerData) = e.fig
figure_backgroundcolor(e::ExplorerData) = e.figure_backgroundcolor
figure_backgroundcolor!(e::ExplorerData, color::RGBA) = e.figure_backgroundcolor = color
window(e::ExplorerData) = e.window

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AXIS 1–4
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ticklabelspace_left(e::ExplorerData) = e.ticklabelspace_left.val
# ticklabelspace_left!(e::ExplorerData, val::Real) = e.ticklabelspace_left[] = val
# ticklabelspace_left_observable(e::ExplorerData) = e.ticklabelspace_left

# ticklabelspace_right(e::ExplorerData) = e.ticklabelspace_right.val
# ticklabelspace_right!(e::ExplorerData, val::Real) = e.ticklabelspace_right[] = val
# ticklabelspace_right_observable(e::ExplorerData) = e.ticklabelspace_right

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AXIS 1
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ax1(e::ExplorerData) = e.ax1

# ax1_height(e::ExplorerData) = e.ax1_height.val

# ax1_font(e::ExplorerData) = e.ax1_font.val
# ax1_font!(e::ExplorerData, font::Any) = e.ax1_font[] = font
# ax1_font_observable(e::ExplorerData) = e.ax1_font

# ax1_fontsize(e::ExplorerData) = e.ax1_fontsize.val
# ax1_fontsize!(e::ExplorerData, fontsize::Int) = e.ax1_fontsize[] = fontsize
# ax1_fontsize_observable(e::ExplorerData) = e.ax1_fontsize

# ax1_textcolor(e::ExplorerData) = e.ax1_textcolor.val
# ax1_textcolor!(e::ExplorerData, color::Symbol) = e.ax1_textcolor[] = color
# ax1_textcolor_observable(e::ExplorerData) = e.ax1_textcolor

# ax1_run_name(e::ExplorerData) = e.ax1_run_name.val
# ax1_run_name!(e::ExplorerData, string::AbstractString) = e.ax1_run_name[] = string
# ax1_run_name_observable(e::ExplorerData) = e.ax1_run_name

# ax1_left_linecolor(e::ExplorerData) = e.ax1_left_linecolor.val
# ax1_left_linecolor!(e::ExplorerData, color::Symbol) = e.ax1_left_linecolor[] = color
# ax1_left_linecolor_observable(e::ExplorerData) = e.ax1_left_linecolor

# ax1_right_linecolor(e::ExplorerData) = e.ax1_right_linecolor.val
# ax1_right_linecolor!(e::ExplorerData, color::Symbol) = e.ax1_right_linecolor[] = color
# ax1_right_linecolor_observable(e::ExplorerData) = e.ax1_right_linecolor

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AXIS 2a
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ax2a(e::ExplorerData) = e.ax2a

# ax2a_height(e::ExplorerData) = e.ax2a_height.val

# ax2a_font(e::ExplorerData) = e.ax2a_font.val
# ax2a_font!(e::ExplorerData, font::Any) = e.ax2a_font[] = font
# ax2a_font_observable(e::ExplorerData) = e.ax2a_font

# ax2a_fontsize(e::ExplorerData) = e.ax2a_fontsize.val
# ax2a_fontsize!(e::ExplorerData, fontsize::Real) = e.ax2a_fontsize[] = fontsize
# ax2a_fontsize_observable(e::ExplorerData) = e.ax2a_fontsize

# ax2a_stored_data_pos_x_y(e::ExplorerData) = e.ax2a_stored_data_pos_x_y.val
# ax2a_stored_data_pos_x_y_observable(e::ExplorerData) = e.ax2a_stored_data_pos_x_y


# ax2a_stored_data_marker(e::ExplorerData) = e.ax2a_stored_data_marker.val
# ax2a_stored_data_marker!(e::ExplorerData, marker::Real) = e.ax2a_stored_data_marker[] = marker
# ax2a_stored_data_marker_observable(e::ExplorerData) = e.ax2a_stored_data_marker

# ax2a_stored_data_markersize(e::ExplorerData) = e.ax2a_stored_data_markersize.val
# ax2a_stored_data_markersize!(e::ExplorerData, size::Real) = e.ax2a_stored_data_markersize[] = size
# ax2a_stored_data_markersize_observable(e::ExplorerData) = e.ax2a_stored_data_markersize

# ax2a_stored_data_color(e::ExplorerData) = e.ax2a_stored_data_color.val
# ax2a_stored_data_color!(e::ExplorerData, color::Real) = e.ax2a_stored_data_color[] = color
# ax2a_stored_data_color_observable(e::ExplorerData) = e.ax2a_stored_data_color

# ax2a_stored_data_plot(e::ExplorerData) = e.ax2a_stored_data_plot

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AXIS 2b - general
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ax2b_xaxis_required(e::ExplorerData) = e.ax2b_xaxis_required.val
# ax2b_xaxis_required!(e::ExplorerData, val::Bool) = e.ax2b_xaxis_required[] = val
# ax2b_xaxis_required_observable(e::ExplorerData) = e.ax2b_xaxis_required

# ax2b_xaxis_datatype(e::ExplorerData) = e.ax2b_xaxis_datatype.val
# ax2b_xaxis_datatype!(e::ExplorerData, datatype::Symbol) = e.ax2b_xaxis_datatype[] = datatype
# ax2b_xaxis_datatype_observable(e::ExplorerData) = e.ax2b_xaxis_datatype_observable

# ax2b_xaxis_minintervalsize(e::ExplorerData) = e.ax2b_xaxis_minintervalsize.val
# ax2b_xaxis_minintervalsize!(e::ExplorerData, minintervalsize::Real) = e.ax2b_xaxis_minintervalsize[] = minintervalsize
# ax2b_xaxis_minintervalsize_observable(e::ExplorerData) = e.ax2b_xaxis_minintervalsize

# ax2b_xaxis_timeunit(e::ExplorerData) = e.ax2b_xaxis_timeunit.val
# ax2b_xaxis_timeunit!(e::ExplorerData, timeunit::Unitful.TimeUnits) = e.ax2b_xaxis_timeunit[] = timeunit
# ax2b_xaxis_timeunit_observable(e::ExplorerData) = e.ax2b_xaxis_timeunit

# ax2b_xlims(e::ExplorerData) = e.ax2b_xlims.val
# ax2b_xlims!(e::ExplorerData, minmax::NTuple{2, Float64}) = e.ax2b_xlims[] = minmax
# ax2b_xlims_observable(e::ExplorerData) = e.ax2b_xlims

# ax2b_xticks_font(e::ExplorerData) = e.ax2b_xticks_font.val
# ax2b_xticks_font!(e::ExplorerData, font::Any) = e.ax2b_xticks_font[] = font
# ax2b_xticks_font_observable(e::ExplorerData) = e.ax2b_xticks_font

# ax2b_xticks_fontsize(e::ExplorerData) = e.ax2b_xticks_fontsize.val
# ax2b_xticks_fontsize!(e::ExplorerData, fontsize::Real) = e.ax2b_xticks_fontsize[] = fontsize
# ax2b_xticks_fontsize_observable(e::ExplorerData) = e.ax2b_xticks_fontsize

# ax2b_yticks_font(e::ExplorerData) = e.ax2b_yticks_font.val
# ax2b_yticks_font!(e::ExplorerData, font::Any) = e.ax2b_yticks_font[] = font
# ax2b_yticks_font_observable(e::ExplorerData) = e.ax2b_yticks_font

# ax2b_yticks_fontsize(e::ExplorerData) = e.ax2b_yticks_fontsize.val
# ax2b_yticks_fontsize!(e::ExplorerData, fontsize::Real) = e.ax2b_yticks_fontsize[] = fontsize
# ax2b_yticks_fontsize_observable(e::ExplorerData) = e.ax2b_yticks_fontsize

# ax2b_pixelwidth(e::ExplorerData) = e.ax2b_pixelwidth.val
# ax2b_pixelwidth_observable(e::ExplorerData) = e.ax2b_pixelwidth

# ax2b_pixelheight(e::ExplorerData) = e.ax2b_pixelheight.val
# ax2b_pixelheight_observable(e::ExplorerData) = e.ax2b_pixelheight

# ax2b_yaxes_extra_stretchfactor(e::ExplorerData) = e.ax2b_yaxes_extra_stretchfactor.val
# ax2b_yaxes_extra_stretchfactor!(e::ExplorerData, stretchfactor::Real) = e.ax2b_yaxes_extra_stretchfactor[] = stretchfactor
# ax2b_yaxes_extra_stretchfactor_observable(e::ExplorerData) = e.ax2b_yaxes_extra_stretchfactor

# ax2b_stretchfactor_changespeed_slow(e::ExplorerData) = e.ax2b_stretchfactor_changespeed_slow
# ax2b_stretchfactor_changespeed_slow!(e::ExplorerData, factor::Real) = e.ax2b_stretchfactor_changespeed_slow = factor

# ax2b_stretchfactor_changespeed_fast(e::ExplorerData) = e.ax2b_stretchfactor_changespeed_fast
# ax2b_stretchfactor_changespeed_fast!(e::ExplorerData, factor::Real) = e.ax2b_stretchfactor_changespeed_fast = factor

# ax2b_linewidth_for_current_runset(e::ExplorerData) = e.ax2b_linewidth_for_current_runset.val
# ax2b_linewidth_for_current_runset!(e::ExplorerData, width::Float64) = e.ax2b_linewidth_for_current_runset[] = width
# ax2b_linewidth_for_current_runset_observable(e::ExplorerData) = e.ax2b_linewidth_for_current_runset

# ax2b_linewidth_for_non_current_runset(e::ExplorerData) = e.ax2b_linewidth_for_non_current_runset.val
# ax2b_linewidth_for_non_current_runset!(e::ExplorerData, width::Float64) = e.ax2b_linewidth_for_non_current_runset[] = width
# ax2b_linewidth_for_non_current_runset_observable(e::ExplorerData) = e.ax2b_linewidth_for_non_current_runset

# ax2b_crosshair_pos_x(e::ExplorerData) = e.ax2b_crosshair_pos_x.val
# ax2b_crosshair_pos_x!(e::ExplorerData, pos::Real) = e.ax2b_crosshair_pos_x[] = pos
# ax2b_crosshair_pos_x_observable(e::ExplorerData) = e.ax2b_crosshair_pos_x

# ax2b_crosshair_pos_y(e::ExplorerData) = e.ax2b_crosshair_pos_y.val
# ax2b_crosshair_pos_y!(e::ExplorerData, pos::Real) = e.ax2b_crosshair_pos_y[] = pos
# ax2b_crosshair_pos_y_observable(e::ExplorerData) = e.ax2b_crosshair_pos_y

# ax2b_crosshair_color(e::ExplorerData) = e.ax2b_crosshair_color.val
# ax2b_crosshair_color!(e::ExplorerData, color::Symbol) = e.ax2b_crosshair_color[] = color
# ax2b_crosshair_color_observable(e::ExplorerData) = e.ax2b_crosshair_color

# ax2b_crosshair_linewidth(e::ExplorerData) = e.ax2b_crosshair_linewidth.val
# ax2b_crosshair_linewidth!(e::ExplorerData, width::Real) = e.ax2b_crosshair_linewidth[] = width
# ax2b_crosshair_linewidth_observable(e::ExplorerData) = e.ax2b_crosshair_linewidth

# ax2b_dragstart(e::ExplorerData) = e.ax2b_dragstart

# ax2b_crosshair_label_font(e::ExplorerData) = e.ax2b_crosshair_label_font.val
# ax2b_crosshair_label_font!(e::ExplorerData, font::Any) = e.ax2b_crosshair_label_font[] = font
# ax2b_crosshair_label_font_observable(e::ExplorerData) = e.ax2b_crosshair_label_font

# ax2b_crosshair_label_fontsize(e::ExplorerData) = e.ax2b_crosshair_label_fontsize.val
# ax2b_crosshair_label_fontsize!(e::ExplorerData, fontsize::Real) = e.ax2b_crosshair_label_fontsize[] = fontsize
# ax2b_crosshair_label_fontsize_observable(e::ExplorerData) = e.ax2b_crosshair_label_fontsize

# ax2b_crosshair_label_positions(e::ExplorerData) = e.ax2b_crosshair_label_positions.val
# ax2b_crosshair_label_positions_observable(e::ExplorerData) = e.ax2b_crosshair_label_positions

# ax2b_crosshair_label_texts(e::ExplorerData) = e.ax2b_crosshair_label_texts.val
# ax2b_crosshair_label_texts_observable(e::ExplorerData) = e.ax2b_crosshair_label_texts

# ax2b_crosshair_label_alignments(e::ExplorerData) = e.ax2b_crosshair_label_alignments.val
# ax2b_crosshair_label_alignments_observable(e::ExplorerData) = e.ax2b_crosshair_label_alignments

# ax2b_crosshair_label_offsets(e::ExplorerData) = e.ax2b_crosshair_label_offsets.val
# ax2b_crosshair_label_offsets_observable(e::ExplorerData) = e.ax2b_crosshair_label_offsets

# ax2b_crosshair_label_visible(e::ExplorerData) = e.ax2b_crosshair_label_visible.val
# ax2b_crosshair_label_visible!(e::ExplorerData, val::Bool) = e.ax2b_crosshair_label_visible[] = val
# ax2b_crosshair_label_visible_observable(e::ExplorerData) = e.ax2b_crosshair_label_visible

# ax2b_gcmsdata_scanindexselect_start(e::ExplorerData) = e.ax2b_gcmsdata_scanindexselect_start.val
# ax2b_gcmsdata_scanindexselect_start!(e::ExplorerData, val::Union{Nothing, Integer}) = e.ax2b_gcmsdata_scanindexselect_start[] = val
# ax2b_gcmsdata_scanindexselect_start_observable(e::ExplorerData) = e.ax2b_gcmsdata_scanindexselect_start

# ax2b_gcmsdata_scanindexselect_stop(e::ExplorerData) = e.ax2b_gcmsdata_scanindexselect_stop.val
# ax2b_gcmsdata_scanindexselect_stop!(e::ExplorerData, val::Union{Nothing, Integer}) = e.ax2b_gcmsdata_scanindexselect_stop[] = val
# ax2b_gcmsdata_scanindexselect_stop_observable(e::ExplorerData) = e.ax2b_gcmsdata_scanindexselect_stop

# ax2b_gcmsdata_scanindexselect_current_range(e::ExplorerData) = e.ax2b_gcmsdata_scanindexselect_current_range.val
# ax2b_gcmsdata_scanindexselect_current_range!(e::ExplorerData, val::Union{Nothing, UnitRange{Int64}}) = e.ax2b_gcmsdata_scanindexselect_current_range[] = val
# ax2b_gcmsdata_scanindexselect_current_range_observable(e::ExplorerData) = e.ax2b_gcmsdata_scanindexselect_current_range

# ax2b_gcmsdata_scanindexselect_mode(e::ExplorerData) = e.ax2b_gcmsdata_scanindexselect_mode.val
# ax2b_gcmsdata_scanindexselect_mode!(e::ExplorerData, val::Union{Nothing, Symbol}) = e.ax2b_gcmsdata_scanindexselect_mode[] = val
# ax2b_gcmsdata_scanindexselect_mode_observable(e::ExplorerData) = e.ax2b_gcmsdata_scanindexselect_mode

# ax2b_gcmsdata_scanindexselect_add(e::ExplorerData) = e.ax2b_gcmsdata_scanindexselect_add.val
# ax2b_gcmsdata_scanindexselect_add_observable(e::ExplorerData) = e.ax2b_gcmsdata_scanindexselect_add

# ax2b_gcmsdata_scanindexselect_subtract(e::ExplorerData) = e.ax2b_gcmsdata_scanindexselect_subtract.val
# ax2b_gcmsdata_scanindexselect_subtract_observable(e::ExplorerData) = e.ax2b_gcmsdata_scanindexselect_subtract

# ax2b_gcmsdata_scantimes_add(e::ExplorerData) = e.ax2b_gcmsdata_scantimes_add.val
# ax2b_gcmsdata_scantimes_add_observable(e::ExplorerData) = e.ax2b_gcmsdata_scantimes_add

# ax2b_gcmsdata_scantimes_add_linewidth(e::ExplorerData) = e.ax2b_gcmsdata_scantimes_add_linewidth.val
# ax2b_gcmsdata_scantimes_add_linewidth!(e::ExplorerData, width::Real) = e.ax2b_gcmsdata_scantimes_add_linewidth[] = width
# ax2b_gcmsdata_scantimes_add_linewidth_observable(e::ExplorerData) = e.ax2b_gcmsdata_scantimes_add_linewidth

# ax2b_gcmsdata_scantimes_add_color(e::ExplorerData) = e.ax2b_gcmsdata_scantimes_add_color.val
# ax2b_gcmsdata_scantimes_add_color!(e::ExplorerData, color::Symbol) = e.ax2b_gcmsdata_scantimes_add_color[] = color
# ax2b_gcmsdata_scantimes_add_color_observable(e::ExplorerData) = e.ax2b_gcmsdata_scantimes_add_color

# ax2b_gcmsdata_scantimes_subtract(e::ExplorerData) = e.ax2b_gcmsdata_scantimes_subtract.val
# ax2b_gcmsdata_scantimes_subtract_observable(e::ExplorerData) = e.ax2b_gcmsdata_scantimes_subtract

# ax2b_gcmsdata_scantimes_subtract_linewidth(e::ExplorerData) = e.ax2b_gcmsdata_scantimes_subtract_linewidth.val
# ax2b_gcmsdata_scantimes_subtract_linewidth!(e::ExplorerData, width::Real) = e.ax2b_gcmsdata_scantimes_subtract_linewidth[] = width
# ax2b_gcmsdata_scantimes_subtract_linewidth_observable(e::ExplorerData) = e.ax2b_gcmsdata_scantimes_subtract_linewidth

# ax2b_gcmsdata_scantimes_subtract_color(e::ExplorerData) = e.ax2b_gcmsdata_scantimes_subtract_color.val
# ax2b_gcmsdata_scantimes_subtract_color!(e::ExplorerData, color::Symbol) = e.ax2b_gcmsdata_scantimes_subtract_color[] = color
# ax2b_gcmsdata_scantimes_subtract_color_observable(e::ExplorerData) = e.ax2b_gcmsdata_scantimes_subtract_color

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AXIS 2b - left y axis
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ax2b_left(e::ExplorerData) = e.ax2b_left

# ax2b_left_required(e::ExplorerData) = e.ax2b_left_required.val
# ax2b_left_required!(e::ExplorerData, val::Bool) = e.ax2b_left_required[] = val
# ax2b_left_required_observable(e::ExplorerData) = e.ax2b_left_required

# ax2b_left_ylims(e::ExplorerData) = e.ax2b_left_ylims.val
# ax2b_left_ylims!(e::ExplorerData, minmax::NTuple{2, Float64}) = e.ax2b_left_ylims[] = minmax
# ax2b_left_ylims_observable(e::ExplorerData) = e.ax2b_left_ylims

# ax2b_left_yaxis_stretchfactor(e::ExplorerData) = e.ax2b_left_yaxis_stretchfactor.val
# ax2b_left_yaxis_stretchfactor!(e::ExplorerData, stretchfactor::Real) = e.ax2b_left_yaxis_stretchfactor[] = stretchfactor
# ax2b_left_yaxis_stretchfactor_observable(e::ExplorerData) = e.ax2b_left_yaxis_stretchfactor

# ax2b_left_prevlimits(e::ExplorerData) = e.ax2b_left_prevlimits

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AXIS 2b - right y axis
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ax2b_right(e::ExplorerData) = e.ax2b_right

# ax2b_right_required(e::ExplorerData) = e.ax2b_right_required.val
# ax2b_right_required!(e::ExplorerData, val::Bool) = e.ax2b_right_required[] = val
# ax2b_right_required_observable(e::ExplorerData) = e.ax2b_right_required

# ax2b_right_ylims(e::ExplorerData) = e.ax2b_right_ylims.val
# ax2b_right_ylims!(e::ExplorerData, minmax::NTuple{2, Float64}) = e.ax2b_right_ylims[] = minmax
# ax2b_right_ylims_observable(e::ExplorerData) = e.ax2b_right_ylims

# ax2b_right_yaxis_stretchfactor(e::ExplorerData) = e.ax2b_right_yaxis_stretchfactor.val
# ax2b_right_yaxis_stretchfactor!(e::ExplorerData, stretchfactor::Real) = e.ax2b_right_yaxis_stretchfactor[] = stretchfactor
# ax2b_right_yaxis_stretchfactor_observable(e::ExplorerData) = e.ax2b_right_yaxis_stretchfactor

# ax2b_right_prevlimits(e::ExplorerData) = e.ax2b_right_prevlimits

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AXIS 3a
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ax3a(e::ExplorerData) = e.ax3a

# ax3a_height(e::ExplorerData) = e.ax3a_height.val

# ax3a_font(e::ExplorerData) = e.ax3a_font.val
# ax3a_font!(e::ExplorerData, font::Any) = e.ax3a_font[] = font
# ax3a_font_observable(e::ExplorerData) = e.ax3a_font

# ax3a_fontsize(e::ExplorerData) = e.ax3a_fontsize.val
# ax3a_fontsize!(e::ExplorerData, fontsize::Real) = e.ax3a_fontsize[] = fontsize
# ax3a_fontsize_observable(e::ExplorerData) = e.ax3a_fontsize

# ax3a_xlims(e::ExplorerData) = e.ax3a_xlims.val
# ax3a_xlims!(e::ExplorerData, minmax::NTuple{2, Float64}) = e.ax3a_xlims[] = minmax
# ax3a_xlims_observable(e::ExplorerData) = e.ax3a_xlims

# ax3a_ylims(e::ExplorerData) = e.ax3a_ylims.val
# ax3a_ylims!(e::ExplorerData, minmax::NTuple{2, Float64}) = e.ax3a_ylims[] = minmax
# ax3a_ylims_observable(e::ExplorerData) = e.ax3a_ylims

# ax3a_extracted_ion_strings(e::ExplorerData) = e.ax3a_extracted_ion_strings.val
# ax3a_extracted_ion_strings_observable(e::ExplorerData) = e.ax3a_extracted_ion_strings

# ax3a_extracted_ion_string_positions(e::ExplorerData) = e.ax3a_extracted_ion_string_positions.val
# ax3a_extracted_ion_string_positions_observable(e::ExplorerData) = e.ax3a_extracted_ion_string_positions

# ax3a_extracted_ion_string_plot(e::ExplorerData) = e.ax3a_extracted_ion_string_plot

# ax3a_ion_intensity_extrema_strings(e::ExplorerData) = e.ax3a_ion_intensity_extrema_strings.val
# ax3a_ion_intensity_extrema_strings_observable(e::ExplorerData) = e.ax3a_ion_intensity_extrema_strings

# ax3a_ion_intensity_extrema_string_positions(e::ExplorerData) = e.ax3a_ion_intensity_extrema_string_positions.val
# ax3a_ion_intensity_extrema_string_positions_observable(e::ExplorerData) = e.ax3a_ion_intensity_extrema_string_positions

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AXIS 3b
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ax3b(e::ExplorerData) = e.ax3b

# ax3b_xlims(e::ExplorerData) = e.ax3b_xlims.val
# ax3b_xlims!(e::ExplorerData, minmax::NTuple{2, Float64}) = e.ax3b_xlims[] = minmax
# ax3b_xlims_observable(e::ExplorerData) = e.ax3b_xlims

# ax3b_ylims(e::ExplorerData) = e.ax3b_ylims.val
# ax3b_ylims!(e::ExplorerData, minmax::NTuple{2, Float64}) = e.ax3b_ylims[] = minmax
# ax3b_ylims_observable(e::ExplorerData) = e.ax3b_ylims

# ax3b_pixelwidth(e::ExplorerData) = e.ax3b_pixelwidth.val
# ax3b_pixelwidth_observable(e::ExplorerData) = e.ax3b_pixelwidth

# ax3b_pixelheight(e::ExplorerData) = e.ax3b_pixelheight.val
# ax3b_pixelheight_observable(e::ExplorerData) = e.ax3b_pixelheight

# ax3b_mode(e::ExplorerData) = e.ax3b_mode.val
# ax3b_mode!(e::ExplorerData, mode::Symbol) = e.ax3b_mode[] = mode
# ax3b_mode_observable(e::ExplorerData) = e.ax3b_mode

# ax3b_ms_ions(e::ExplorerData) = e.ax3b_ms_ions.val

# ax3b_ms_ions_observable(e::ExplorerData) = e.ax3b_ms_ions
# ax3b_ms_lineheights(e::ExplorerData) = e.ax3b_ms_lineheights.val
# ax3b_ms_lineheights_observable(e::ExplorerData) = e.ax3b_ms_lineheights

# ax3b_ms_linecolor(e::ExplorerData) = e.ax3b_ms_linecolor.val
# ax3b_ms_linecolor!(e::ExplorerData, color::Symbol) = e.ax3b_ms_linecolor[] = color
# ax3b_ms_linecolor_observable(e::ExplorerData) = e.ax3b_ms_linecolor

# ax3b_ms_label_positions(e::ExplorerData) = e.ax3b_ms_label_positions.val
# ax3b_ms_label_positions_observable(e::ExplorerData) = e.ax3b_ms_label_positions

# ax3b_ms_label_texts(e::ExplorerData) = e.ax3b_ms_label_texts.val
# ax3b_ms_label_texts_observable(e::ExplorerData) = e.ax3b_ms_label_texts

# ax3b_ion_chromatogram_scantimes(e::ExplorerData) = e.ax3b_ion_chromatogram_scantimes.val
# ax3b_ion_chromatogram_scantimes_observable(e::ExplorerData) = e.ax3b_ion_chromatogram_scantimes

# ax3b_ion_chromatogram_intensities(e::ExplorerData) = e.ax3b_ion_chromatogram_intensities.val
# ax3b_ion_chromatogram_intensities_observable(e::ExplorerData) = e.ax3b_ion_chromatogram_intensities

# ax3b_ion_chromatogram_linewidth(e::ExplorerData) = e.ax3b_ion_chromatogram_linewidth.val
# ax3b_ion_chromatogram_linewidth!(e::ExplorerData, width::Real) = e.ax3b_ion_chromatogram_linewidth[] = width
# ax3b_ion_chromatogram_linewidth_observable(e::ExplorerData) = e.ax3b_ion_chromatogram_linewidth

# ax3b_ion_chromatogram_color(e::ExplorerData) = e.ax3b_ion_chromatogram_color.val
# ax3b_ion_chromatogram_color!(e::ExplorerData, color::Symbol) = e.ax3b_ion_chromatogram_color[] = color
# ax3b_ion_chromatogram_color_observable(e::ExplorerData) = e.ax3b_ion_chromatogram_color

# ax3b_crosshair_pos_x(e::ExplorerData) = e.ax3b_crosshair_pos_x.val
# ax3b_crosshair_pos_x_observable(e::ExplorerData) = e.ax3b_crosshair_pos_x

# ax3b_crosshair_color(e::ExplorerData) = e.ax3b_crosshair_color.val
# ax3b_crosshair_color!(e::ExplorerData, color::Symbol) = e.ax3b_crosshair_color[] = color
# ax3b_crosshair_color_observable(e::ExplorerData) = e.ax3b_crosshair_color.val

# ax3b_crosshair_linewidth(e::ExplorerData) = e.ax3b_crosshair_linewidth.val
# ax3b_crosshair_linewidth!(e::ExplorerData, width::Real) = e.ax3b_crosshair_linewidth[] = width
# ax3b_crosshair_linewidth_observable(e::ExplorerData) = e.ax3b_crosshair_linewidth.val

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AXIS 4a
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ax4a(e::ExplorerData) = e.ax4a

# ax4a_height(e::ExplorerData) = e.ax4a_height.val

# ax4a_font(e::ExplorerData) = e.ax4a_font.val

# ax4a_fontsize(e::ExplorerData) = e.ax4a_fontsize.val
# ax4a_fontsize!(e::ExplorerData, fontsize::Real) = e.ax4a_fontsize[] = fontsize

# ax4a_xlims(e::ExplorerData) = e.ax4a_xlims.val
# ax4a_xlims!(e::ExplorerData, minmax::NTuple{2, Float64}) = e.ax4a_xlims[] = minmax
# ax4a_xlims_observable(e::ExplorerData) = e.ax4a_xlims

# ax4a_ylims(e::ExplorerData) = e.ax4a_ylims.val
# ax4a_ylims!(e::ExplorerData, minmax::NTuple{2, Float64}) = e.ax4a_ylims[] = minmax
# ax4a_ylims_observable(e::ExplorerData) = e.ax4a_ylims

# ax4a_ion_positions(e::ExplorerData) = e.ax4a_ion_positions.val
# ax4a_ion_positions_observable(e::ExplorerData) = e.ax4a_ion_positions

# ax4a_ion_strings(e::ExplorerData) = e.ax4a_ion_strings.val
# ax4a_ion_strings_observable(e::ExplorerData) = e.ax4a_ion_strings

# ax4a_ion_colors(e::ExplorerData) = e.ax4a_ion_colors.val
# ax4a_ion_colors_observable(e::ExplorerData) = e.ax4a_ion_colors

# ax4a_ion_offsets(e::ExplorerData) = e.ax4a_ion_offsets.val
# ax4a_ion_offsets_observable(e::ExplorerData) = e.ax4a_ion_offsets

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AXIS 4b
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ax4b(e::ExplorerData) = e.ax4b

# ax4b_xlims(e::ExplorerData) = e.ax4b_xlims.val
# ax4b_xlims!(e::ExplorerData, minmax::NTuple{2, Float64}) = e.ax4b_xlims[] = minmax
# ax4b_xlims_observable(e::ExplorerData) = e.ax4b_xlims

# ax4b_ylims(e::ExplorerData) = e.ax4b_ylims.val
# ax4b_ylims!(e::ExplorerData, minmax::NTuple{2, Float64}) = e.ax4b_ylims[] = minmax
# ax4b_ylims_observable(e::ExplorerData) = e.ax4b_ylims

# ax4b_xticks_font(e::ExplorerData) = e.ax4b_xticks_font.val
# ax4b_xticks_font!(e::ExplorerData, font::Any) = e.ax4b_xticks_font[] = font
# ax4b_xticks_font_observable(e::ExplorerData) = e.ax4b_xticks_font

# ax4b_xticks_fontsize(e::ExplorerData) = e.ax4b_xticks_fontsize.val
# ax4b_xticks_fontsize!(e::ExplorerData, fontsize::Real) = e.ax4b_xticks_fontsize[] = fontsize
# ax4b_xticks_fontsize_observable(e::ExplorerData) = e.ax4b_xticks_fontsize

# ax4b_pixelwidth(e::ExplorerData) = e.ax4b_pixelwidth.val
# ax4b_pixelwidth_observable(e::ExplorerData) = e.ax4b_pixelwidth

# ax4b_pixelheight(e::ExplorerData) = e.ax4b_pixelheight.val
# ax4b_pixelheight_observable(e::ExplorerData) = e.ax4b_pixelheight

# ax4b_yaxes_extra_stretchfactor(e::ExplorerData) = e.ax4b_yaxes_extra_stretchfactor.val
# ax4b_yaxes_extra_stretchfactor!(e::ExplorerData, factor::Real) = e.ax4b_yaxes_extra_stretchfactor[] = factor
# ax4b_yaxes_extra_stretchfactor_observable(e::ExplorerData) = e.ax4b_yaxes_extra_stretchfactor
# ax4b_stretchfactor_changespeed_slow(e::ExplorerData) = e.ax4b_stretchfactor_changespeed_slow
# ax4b_stretchfactor_changespeed_slow!(e::ExplorerData, factor::Real) = e.ax4b_stretchfactor_changespeed_slow = factor
# ax4b_stretchfactor_changespeed_fast(e::ExplorerData) = e.ax4b_stretchfactor_changespeed_fast
# ax4b_stretchfactor_changespeed_fast!(e::ExplorerData, factor::Real) = e.ax4b_stretchfactor_changespeed_fast = factor

# ax4b_ms_ions(e::ExplorerData) = e.ax4b_ms_ions.val
# ax4b_ms_ions_observable(e::ExplorerData) = e.ax4b_ms_ions

# ax4b_ms_lineheights(e::ExplorerData) = e.ax4b_ms_lineheights.val
# ax4b_ms_lineheights_observable(e::ExplorerData) = e.ax4b_ms_lineheights

# ax4b_ms_linewidth(e::ExplorerData) = e.ax4b_ms_linewidth.val
# ax4b_ms_linewidth!(e::ExplorerData, width::Real) = e.ax4b_ms_linewidth[] = width
# ax4b_ms_linewidth_observable(e::ExplorerData) = e.ax4b_ms_linewidth

# ax4b_ms_linecolor(e::ExplorerData) = e.ax4b_ms_linecolor.val
# ax4b_ms_linecolor!(e::ExplorerData, color::Symbol) = e.ax4b_ms_linecolor[] = color
# ax4b_ms_linecolor_observable(e::ExplorerData) = e.ax4b_ms_linecolor

# ax4b_ms_label_positions(e::ExplorerData) = e.ax4b_ms_label_positions.val
# ax4b_ms_label_positions_observable(e::ExplorerData) = e.ax4b_ms_label_positions

# ax4b_ms_label_texts(e::ExplorerData) = e.ax4b_ms_label_texts.val
# ax4b_ms_label_texts_observable(e::ExplorerData) = e.ax4b_ms_label_texts

# ax4b_ms_label_font(e::ExplorerData) = e.ax4b_ms_label_font.val
# ax4b_ms_label_font!(e::ExplorerData, font::Any) = e.ax4b_ms_label_font[] = font
# ax4b_ms_label_font_observable(e::ExplorerData) = e.ax4b_ms_label_font

# ax4b_ms_label_fontsize(e::ExplorerData) = e.ax4b_ms_label_fontsize.val
# ax4b_ms_label_fontsize!(e::ExplorerData, fontsize::Real) = e.ax4b_ms_label_fontsize[] = fontsize
# ax4b_ms_label_fontsize_observable(e::ExplorerData) = e.ax4b_ms_label_fontsize

# ax4b_ms_label_offset(e::ExplorerData) = e.ax4b_ms_label_offset.val
# ax4b_ms_label_offset!(e::ExplorerData, offset::NTuple{2, Int}) = e.ax4b_ms_label_offset[] = offset
# ax4b_ms_label_offset_observable(e::ExplorerData) = e.ax4b_ms_label_offset

# ax4b_slider_position(e::ExplorerData) = e.ax4b_slider_position.val
# ax4b_slider_position_observable(e::ExplorerData) = e.ax4b_slider_position

# ax4b_slider_linewidth(e::ExplorerData) = e.ax4b_slider_linewidth.val
# ax4b_slider_linewidth!(e::ExplorerData, width::Real) = e.ax4b_slider_linewidth[] = width
# ax4b_slider_linewidth_observable(e::ExplorerData) = e.ax4b_slider_linewidth

# ax4b_slider_linecolor(e::ExplorerData) = e.ax4b_slider_linecolor.val
# ax4b_slider_linecolor!(e::ExplorerData, color::Symbol) = e.ax4b_slider_linecolor[] = color
# ax4b_slider_linecolor_observable(e::ExplorerData) = e.ax4b_slider_linecolor

# ax4b_selected_ion_index(e::ExplorerData) = e.ax4b_selected_ion_index.val
# ax4b_selected_ion_index!(e::ExplorerData, index::Union{Nothing, Integer}) = e.ax4b_selected_ion_index[] = index
# ax4b_selected_ion_index_observable(e::ExplorerData) = e.ax4b_selected_ion_index

# ax4b_selected_ions(e::ExplorerData) = e.ax4b_selected_ions.val
# ax4b_selected_ions_observable(e::ExplorerData) = e.ax4b_selected_ions

# ax4b_selected_ion_positions(e::ExplorerData) = e.ax4b_selected_ion_positionss.val
# ax4b_selected_ion_positions_observable(e::ExplorerData) = e.ax4b_selected_ion_positions

# ax4b_selected_ion_colors(e::ExplorerData) = e.ax4b_selected_ion_colors.val
# ax4b_selected_ion_colors_observable(e::ExplorerData) = e.ax4b_selected_ion_colors

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Fields (and utility functions) related to the storage and selection of runsets.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# runsets(e::ExplorerData) = e.runsets.val
# runsets_observable(e::ExplorerData) = e.runsets

# runsets_count(e::ExplorerData) = length(e.runsets.val)

# runset_index(e::ExplorerData) = e.runset_index.val
# runset_index!(e::ExplorerData, index::Union{Integer, Nothing}) = e.runset_index[] = index
# runset_index_observable(e::ExplorerData) = e.runset_index

# current_runset(e::ExplorerData) = runsets_count(e) > 0 ? runsets(e)[runset_index(e)] : nothing

# current_runset_cleanup(e::ExplorerData) = e.current_runset_cleanup.val
# current_runset_cleanup!(e::ExplorerData, index::Union{Union, Integer}) = e.current_runset_cleanup[] = index
# current_runset_cleanup_observable(e::ExplorerData) = e.current_runset_cleanup


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Fields related to the import of runsets
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# filereader(e::ExplorerData) = e.filereader
# filereader!(e::ExplorerData, reader::Union{Function, Nothing}) = e.filereader = reader
