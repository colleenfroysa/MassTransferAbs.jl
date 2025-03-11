mutable struct HenrysPoly{T}
    C₁::T
    C₂::T
    C₃::T
    C₄::T
end

mutable struct PackingMaterial{T} # Ralu pak 
    a::T
    ϵ::T 
    C_L::T
    C_V::T
end

mutable struct SpecificInterfaceArea{T} # aₚₕ
    A₁::T
    A₂::T
    A₃::T
    A₄::T
    A₅::T
end