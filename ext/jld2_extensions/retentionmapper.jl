struct RetentionMapperJLD2V1{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14,T15, T16}
    rA::T1
    rA_unit::Union{Nothing, T2}
    rA_min::T3
    rA_max::T4
    rA_norm_min::T5
    rA_norm_max::T6
    rB::T7
    rB_unit::Union{Nothing, T8}
    rB_min::T9
    rB_max::T10
    rB_pred_min::T11
    rB_pred_max::T12
    rB_pred_norm_min::T13
    rB_pred_norm_max::T14
    knots::T15
    coefs::T16
    extras::Dict{String, Any}
end

JLD2.writeas(::Type{JuChrom.RetentionMapper{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17}}) where {T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17} =
    RetentionMapperJLD2V1{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16}

function JLD2.wconvert(
    ::Type{RetentionMapperJLD2V1{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, 
        T16}},
    rm::JuChrom.RetentionMapper{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, 
        T16, T17}
    ) where {T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14,T15, T16, T17}

    RetentionMapperJLD2V1{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16}(
        rm.rA,
        rm.rA_unit,
        rm.rA_min,
        rm.rA_max,
        rm.rA_norm_min,
        rm.rA_norm_max,
        rm.rB,
        rm.rB_unit,
        rm.rB_min,
        rm.rB_max,
        rm.rB_pred_min,
        rm.rB_pred_max,
        rm.rB_pred_norm_min,
        rm.rB_pred_norm_max,
        rm.knots,
        rm.coefs,
        rm.extras
    )
end

function JLD2.rconvert(
    ::Type{JuChrom.RetentionMapper{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, 
        T16, T17}}, 
    rm::RetentionMapperJLD2V1{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16}
    ) where {T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17}
    
    JuChrom.RetentionMapper(
        rm.rA,
        rm.rA_unit,
        rm.rA_min,
        rm.rA_max,
        rm.rA_norm_min,
        rm.rA_norm_max,
        rm.rB,
        rm.rB_unit,
        rm.rB_min,
        rm.rB_max,
        rm.rB_pred_min,
        rm.rB_pred_max,
        rm.rB_pred_norm_min,
        rm.rB_pred_norm_max,
        rm.knots,
        rm.coefs,
        BSplineKit.Spline(BSplineKit.RecombinedBSplineBasis(BSplineKit.BSplineBasis(
            BSplineKit.BSplineOrder(4), rm.knots), BSplineKit.Natural()), rm.coefs),
        rm.extras
    )
end
