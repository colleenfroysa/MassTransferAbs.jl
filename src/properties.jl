using UnPack

function flux(packing_info::A, K::T, ΔC::T, correction_factor::T) where {T <: Number, A <: PackingMaterial}
    Aw = Aw(packing_info, correction_factor)  
    return K * Aw * ΔC
end


function overall_mass_transfer(liquid_phase::F, vapour_phase::F, packing::A, t::T, x::T, y::T, u_L::T, u_v::T, v_L::T, v_v::T) where {T <: Number, A <: PackingMaterial, F <: SimpleFluidMedium} 
    Kₗ = mass_transfer(liquid_phase, packing, t, x, u_L, v_L)
    Kᵥ = mass_transfer(vapour_phase, packing, t, y, u_v, v_v)
    E = enhancement_factor()
end

function mass_transfer(phase::F, packing::A, t::T, x::T, u::T, v::T) where {T <: Number, A <: PackingMaterial, F <: SimpleFluidMedium} 
    @unpack a, ϵ, C_L, C_V = packing
    D = diffusivity(phase, t, x)  # Calls the right diffusivity function automatically
    h_L = liquid_holdup(packing, phase, u)
    d_h = hydraulic_diameter(packing)
end

function liquid_holdup(packing::A, phase::F, u_L::T) where {T <: Number, A <: PackingMaterial, F <: SimpleFluidMedium}
    @unpack a, ϵ, C_L, C_V = packing
    η_L = η_unloaded(phase.data.viscosity_data)  # Get viscosity from `SimpleFluidMedium`
    ρ_L = ρ_unloaded(phase.data.density_data) # Get density from `SimpleFluidMedium`
    g = 9.81  # Gravity
    return (12 * (1 / g) * (η_L / ρ_L) * u_L * a^2)^(1/3)
end

function hydraulic_diameter(packing::A) where A <: PackingMaterial
    @unpack a, ϵ, C_L, C_V = packing
    return 4 * (ϵ / a)
end

function specific_interface_area(packing::A, specific_area::B, d_h::T, Re_L::T, We_L::T, Fr_L::T) where {T <: Number, A <: PackingMaterial, B <: SpecificInterfaceArea}
    @unpack a, ϵ, C_L, C_V = packing
    @unpack A₁, A₂, A₃, A₄, A₅ = specific_area
    return A₁ * (a * d_h)^A₂ * Re_L^A₃ * We_L^A₄ * Fr_L^A₅ * a
end

function dimensionless_numbers(u::T, ρ::T, d_h::T, η::T, σ::T) where T <: Number
    g = 9.81  # Gravity
    Re = (u * d_h) / η    # Reynolds number
    We = (u^2 * ρ * d_h) / σ  # Weber number
    Fr = u^2 / (g * d_h)   # Froude number
    return Re, We, Fr
end

function enhancement_factor()
    k_CO₂ = kinetic_rate_constant(x_MEA, t, Θ, α)
    Dₗ = diffusivity(phase, t, x)
    
end

function kinetic_rate_constant(x_MEA::T, t::T, Θ::T, α::T) where T <: Number
    k_2_star = k₂_star(Val(x_MEA), t)
    return k_2_star * (1 / (exp(1 / (Θ - α))))
end

function k₂_star(::Val{1}, t::T) where T <: Number
    return 8.87e5 * exp(-3458 / t)
end

function k₂_star(::Val{5}, t::T) where T <: Number
    return 4.396e6 * exp(-3693 / t)
end

function diffusivity(liquid_phase::F, t::T, x::T) where {T <: Number, F <: SimpleFluidMedium} 
    # Compute diffusivity for liquid
    return 2.0398 * 10^(-9)  # Example value 
end

function diffusivity(vapour_phase::F, t::T, x::T) where {T <: Number, F <: SimpleFluidMedium} 
    # Compute diffusivity for vapour
    return 1.933 * 10^(-5)  # Example value
end

# -------------------------------------

function Aw(packing_info::A, correction_factor::T) where {T <: Number, A <: PackingMaterial}
    @unpack a, ϵ, C_L, C_V = packing_info
    return a * correction_factor
end

function henrys_constant(henrys_constant_model::HenrysPoly, t::T) where T <: Number
    @unpack C₁, C₂, C₃, C₄ = henrys_constant_model
    return C₁ + C₂ / t + C₃ / ln(t) + C₄ / t 
end

function v_V(μ_V::T, ρ_V::T) where T <: Number
    return μ_V / ρ_V
end