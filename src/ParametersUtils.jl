module ParametersUtils
using TOML

Base.@kwdef struct Parameters
    TC::Float64
    Gamma::Float64
    gamma_spr::Float64
    gamma_r::Float64
    c::Float64
    reduced_planck::Float64
    T::Float64
    omega_hfs::Float64
    omega_r1::Float64
    Delta_omega::Float64
    omega_r2::Float64
    a_m::Float64
    a_Tem::Float64
    m::Float64
    V_T::Float64
    V_Tr::Float64
    Gamma_L::Float64
    k1::Float64
    k2::Float64
    d31::Float64
    d32::Float64
    na::Float64
    Omega1::Float64
    Omega2::Float64
    V1_0::ComplexF64
    V2_0::ComplexF64
    a::Float64
end

function load_from_file(file)
    data = TOML.parse(file)
    laser = data["laser"]
    constants = data["constants"]

    TC = laser["TC"]
    Gamma = laser["Gamma"]
    gamma_spr = laser["gamma_spr"]
    gamma_r = laser["gamma_r"]


    c = constants["c"]
    reduced_planck = constants["reduced_planck"]
    omega_hfs = constants["omega_hfs"]
    omega_r1 = constants["omega_r1"]
    a_m = constants["a_m"]
    a_Tem = constants["a_Tem"]
    Omega1 = constants["Omega1"]
    Omega2 = constants["Omega2"]
    V1_0 = Complex(constants["V1_0"]...)
    V2_0 = Complex(constants["V2_0"]...)
    a = constants["a"]

    T = (TC + 273.0) / 11604.5
    Delta_omega = omega_hfs
    omega_r2 = omega_r1 - Delta_omega
    m = a_m * a_Tem
    V_T = sqrt(2 * T / m)
    V_Tr = V_T * c
    Gamma_L = Gamma / 2.0
    k1 = omega_r1 / c
    k2 = omega_r2 / c
    d31 = sqrt(3 * gamma_spr * reduced_planck * c^3 / (2 * omega_r1^3))
    d32 = sqrt(3 * gamma_spr * reduced_planck * c^3 / (2 * omega_r2^3))
    na = exp(-0.0003 * TC^2 + 0.119 * TC + 20.4)

    Parameters(
        TC=TC,
        Gamma=Gamma,
        gamma_spr=gamma_spr,
        gamma_r=gamma_r,
        c=c,
        reduced_planck=reduced_planck,
        T=T,
        omega_hfs=omega_hfs,
        omega_r1=omega_r1,
        Delta_omega=Delta_omega,
        omega_r2=omega_r2,
        a_m=a_m,
        a_Tem=a_Tem,
        m=m,
        V_T=V_T,
        V_Tr=V_Tr,
        Gamma_L=Gamma_L,
        k1=k1,
        k2=k2,
        d31=d31,
        d32=d32,
        na=na,
        Omega1=Omega1,
        Omega2=Omega2,
        V1_0=V1_0,
        V2_0=V2_0,
        a=a,
    )
end

function dump_to_file(file, parameters::Parameters)
    println(
        file,
        """
[laser]
TC = $(parameters.TC)
Gamma = $(parameters.Gamma)
gamma_spr = $(parameters.gamma_spr)
gamma_r = $(parameters.gamma_r)
    
[constants]
c = $(parameters.c)
reduced_planck = $(parameters.reduced_planck)
omega_hfs = $(parameters.omega_hfs)
omega_r1 = $(parameters.omega_r1)
a_m = $(parameters.a_m)
a_Tem = $(parameters.a_Tem)
Omega1 = $(parameters.Omega1)
Omega2 = $(parameters.Omega2)
V1_0 = [$(parameters.V1_0.re), $(parameters.V1_0.im)]
V2_0 = [$(parameters.V2_0.re), $(parameters.V2_0.im)]
a = $(parameters.a)

[computed]
T = $(parameters.T)
Delta_omega = $(parameters.Delta_omega)
omega_r2 = $(parameters.omega_r2)
m = $(parameters.m)
V_T = $(parameters.V_T)
V_Tr = $(parameters.V_Tr)
Gamma_L = $(parameters.Gamma_L)
k1 = $(parameters.k1)
k2 = $(parameters.k2)
d31 = $(parameters.d31)
d32 = $(parameters.d32)
na = $(parameters.na)
"""
    )
end

end
