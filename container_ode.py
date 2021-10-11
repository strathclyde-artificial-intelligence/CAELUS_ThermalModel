#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 19:13:36 2021

@author: gianlucafilippi
"""

module Container_2

export CalculatePayloadProp, CalculatePCMProp, CalculateContainerProp, CalculateAirProp
export CalculateMatrices
export container_thermalnetwork!

using LinearAlgebra

# ------ General Structs -------

"""
    Container

"""
struct ContainerProp{TF}
    x_in::TF
    y_in::TF
    z_in::TF
    x_out::TF
    y_out::TF
    z_out::TF
    k::TF      # conductivity
    ϵ::TF      # emissivity
    cp::TF     # specific heat
    ρ::TF      # density
    Ain_cond::TF
    Ain_conv::TF
    Aout::TF
    V::TF
end


"""
    Payload

"""
struct PayloadProp{TF}
    x::TF
    y::TF
    z::TF
    k::TF      # conductivity
    ϵ::TF      # emissivity
    cp::TF     # specific heat
    ρ::TF      # density
    A_cond::TF
    A_conv::TF
    V::TF
end



"""
    Phase Change Material

"""
struct PCMProp{TF}
    x::TF
    y::TF
    z::TF
    k_s::TF      # conductivity
    k_l::TF      # conductivity
    ϵ::TF      # emissivity
    cp_s::TF     # specific heat
    cp_l::TF     # specific heat
    ρ::TF      # density
    Tm1::TF
    Tm2::TF
    Hlatent::TF
    A_cond::TF
    A_conv::TF
    V::TF
end


"""
    Internal air

"""
struct AirProp{TF}
    cp::TF     # specific heat
    ρ::TF      # density
    h_NC::TF   # natural convection coefficient
    h_FC::TF   # forced convection coefficient
    V::TF
end




# ------------- helper functions --------------


"""
    CalculatePayloadProp(state)

    x[1]: length x [m]
    x[2]: width y [m]
    x[3]: high z [m]
    x[4]: conductivity k [W m-1 K-1]
    x[5]: emissivity ϵ
    x[6]: cp [J kg-1 K-1]
    x[7]: density ρ [kg m-3]
"""
function CalculatePayloadProp(x)

    A_conv = 2* ( x[1]*x[3] + x[2]*x[3] )
    A_cond = x[1] * x[2]
    V = prod(x[1:3])

    x_tot = [x; A_cond; A_conv; V]
    out =  PayloadProp(x_tot...)
end







"""
    CalculatePCMProp(state)
sfab fdbabe
"""
function CalculatePCMProp(x,pay)

    A_cond = pay.A_cond
    A_conv = 2* ( x[1]*x[2] + x[1]*x[3] + x[2]*x[3] ) - pay.A_cond
    V = prod(x[1:3])

    x_tot = [x; A_cond; A_conv; V]
    out =  PCMProp(x_tot...)
end


"""
    CalculateContainerProp(state)
sfab fdbabe
"""
function CalculateContainerProp(x,pay)

    A_in = 2 * ( x[1]*x[2] + x[1]*x[3] + x[2]*x[3] )
    A_out = 2 * ( x[4]*x[5] + x[4]*x[6] + x[5]*x[6] )
    A_in_cond = pay.A_cond
    A_in_conv = A_in - A_in_cond
    V = prod(x[4:6]) - prod(x[1:3])

    x_tot = [x; A_in_cond; A_in_conv; A_out; V]
    out =  ContainerProp(x_tot...)
end



"""
    CalculateAirProp(state)
sfab fdbabe
"""
function CalculateAirProp(x, con, pay, pcm)

    V = con.V - pay.V -  pcm.V

    x_tot = [x; V]
    out =  AirProp(x_tot...)
end





"""
    CalculateAreas(state)
sfab fdbabe
"""
function CalculateAreas(con_prop, pay_prop, pcm_prop)


    L_brick = [pay_prop.x, pay_prop.y, pay_prop.z]
    L_PCM  = [pcm_prop.x, pcm_prop.y, pcm_prop.z]
    L_in_container = [con_prop.x_in, con_prop.y_in, con_prop.z_in]
    L_out_container = [con_prop.x_out, con_prop.y_out, con_prop.z_out]



    A_brick = 2*(L_brick[1] * L_brick[2]) + 2*(L_brick[1] * L_brick[3]) + 2*(L_brick[2] * L_brick[3])              # total external area container
    A_brick_conv = 2* ( L_brick[1]*L_brick[3] + L_brick[2]*L_brick[3] )
    A_brick_cond = (L_brick[1] * L_brick[2])


    A_PCM = 2*(L_PCM[1] * L_PCM[2]) + 2*(L_PCM[1] * L_PCM[3]) + 2*(L_PCM[2] * L_PCM[3])
    A_PCM_conv = A_PCM - A_brick_cond
    A_PCM_cond = (L_PCM[1] * L_PCM[2])


    A_container_out = 2*(L_out_container[1] * L_out_container[2]) + 2*(L_out_container[1] * L_out_container[3]) + 2*(L_out_container[2] * L_out_container[3])
    A_container_in  = 2*(L_in_container[1] * L_in_container[2]) + 2*(L_in_container[1] * L_out_container[3]) + 2*(L_in_container[2] * L_out_container[3])
    # A_container_in_convection = A_container_in - Lx_brick*Ly_brick


    A_container_in_cond = L_brick[1] * L_brick[2]
    A_container_mean = (A_container_out+A_container_in) / 2
    A_container_in_conv = A_container_in - A_brick_cond

    return A_brick_cond, A_brick_conv, A_PCM_cond, A_PCM_conv, A_container_in_cond, A_container_in_conv, A_container_out

end



"""
    CalculateVolumes(state)
sfab fdbabe
"""
function CalculateVolumes(con_prop, pay_prop, pcm_prop)


    L_brick = [pay_prop.x, pay_prop.y, pay_prop.z]
    L_PCM  = [pcm_prop.x, pcm_prop.y, pcm_prop.z]
    L_in_container = [con_prop.x_in, con_prop.y_in, con_prop.z_in]
    L_out_container = [con_prop.x_out, con_prop.y_out, con_prop.z_out]



    V_brick = prod( L_brick )
    V_PCM = prod( L_PCM )

    V_out_container = prod( L_out_container )
    V_in_container = prod( L_in_container )
    V_Dx_container = V_out_container - V_in_container

    V_air = V_in_container - V_brick -  V_PCM

    return V_brick, V_PCM, V_Dx_container, V_air

end




"""
    SpreadingResistance(state)

"""
function SpreadingResistance(r1, r2, t, k, h)

    ϵ = r1/r2
    τ = t / r2
    Bi = (h * r2) / k
    λ = pi + 1/(1/ (ϵ*sqrt(pi)) )
    ϕ = (tanh(λ*τ) + λ/Bi) / (1+ λ/Bi*tanh(λ*τ))
    ψ = (ϵ*τ)/sqrt(pi) + 1/sqrt(pi)*(1-ϵ)*ϕ

    Rsp = ψ / (k*r1*sqrt(pi))
    return Rsp
end



"""
    CalculateMatrices(state)
sfab fdbabe
"""
function CalculateMatrices(con_prop, pay_prop, pcm_prop, air_prop)


    # ------------- initialisation ------------------------
    # matrix initialisation
    H = zeros(6,6)  # convection
    K = zeros(6,6)  # conduction
    R = zeros(6,6)  # radiation
    Q = zeros(6)    # external heat source

    # M = zeros(5)    # vector of masses
    C = zeros(6)    # vector of heat capacity
    # -----------------------------------------------------




    # --------------- geometry calculationns ----------------------------
    A_payload_cond, A_payload_conv, A_PCM_cond, A_PCM_conv, A_container_in_cond, A_container_in_conv, A_container_out = CalculateAreas(con_prop, pay_prop, pcm_prop)
    V_payload, V_PCM, V_container, V_air = CalculateVolumes(con_prop, pay_prop, pcm_prop)
    # -------------------------------------------------------



    # -------------- heat capcity ------------------------
    # node 1 : external wall / payload
    C[1] = V_container * con_prop.ρ * con_prop.cp   # J K-1

    # node 2 : internal wall / payload
    C[2] = V_container * con_prop.ρ * con_prop.cp   # J K-1

    # node 3 : payload
    C[3] = V_payload * pay_prop.ρ * pay_prop.cp

    # node 4 : PCM

    # node 5 : air
    C[5] = V_air * air_prop.ρ * air_prop.cp
    # -------------------------------------------------

    Δx_container = ( con_prop.z_out - con_prop.z_in ) / 2



    # -------- spreading resistances--------

    # from "square" to "circular" geometry
    r_contianer = sqrt(A_container_out / pi)
    r_brick = sqrt(A_payload_cond / pi)
    r_PCM = sqrt(A_PCM_cond / pi)

    R_c2b = SpreadingResistance(r_brick, r_contianer, Δx_container, con_prop.k, air_prop.h_NC)
    R_PCM2b = SpreadingResistance(r_brick, r_PCM, pcm_prop.z, pcm_prop.k_s, air_prop.h_NC)

    # R_c2b=0
    # R_PCM2b=0
    # ------------------------------------------


    # ----------------- conductivity ------------------------

    # node 1 (external wall) <-> node 2 (internal wall/payload)
    K[1,2] = con_prop.k * A_container_out / Δx_container
    K[2,1] = K[1,2]

    # node 2 (internal wall/payload)  <-> node 3 (payload)
    R_23 = ( pay_prop.k * A_payload_cond / (pay_prop.z/2) )^(-1)
    R_tot = R_23+R_c2b
    K_tot = R_tot^-1
    K[3,2] = K_tot
    K[2,3] = K[3,2]#R_23^-1

    # node 3 (payload)  <-> node 4 (PCM)                 N.B.! ONLY PAYLOAD PART
    R_34 = ( pay_prop.k * A_PCM_cond / (pay_prop.z/2) )^-1
    R_tot = R_34+R_PCM2b
    K[4,3] = R_34^-1
    K[3,4] = K[4,3]#R_tot^(-1)
    Kpartial = K[4,3]
    # ---------------------------------------------------------



    # --------- convection -------------
    # node 1 (internal wall)  <-> node 6 (external air)
    H[1,6] = air_prop.h_NC * A_container_out
    H[6,1] = H[1,6]

    # node 2 (internal wall)  <-> node 5 (air)
    H[2,5] = air_prop.h_NC * A_container_in_conv
    H[5,2] = H[2,5]

    # node 3 (payload)  <-> node 5 (air)
    H[3,5] = air_prop.h_NC * A_payload_conv
    H[5,3] = H[3,5]

    # node 4 (PCM)  <-> node 5 (air)
    H[4,5] = air_prop.h_NC * A_PCM_conv
    H[5,4] = H[4,5]
    Hpartial = H[4,5]
    # ----------------------------------

    return C, K, H, R, Q, Kpartial, Hpartial
end

# ------------------------------------------------------------------------------




# --------------------- main function (public) -------------------

function container_thermalnetwork!(dT, T, params, t)

    C, K, H, R, Q, Kpartial, Hpartial, pcm_prop, atmmodel = params

    # ambient temperature from Atmosferic model
    T[6] = atmmodel(t)


    # PCM interpolation
    K4 = []
    L = (T[4]-pcm_prop.Tm1) / (pcm_prop.Tm2-pcm_prop.Tm1)
    # PCM properties
    if L < 0
        Ln = 0
        Mass_pcm = pcm_prop.ρ * pcm_prop.V
        cp = pcm_prop.cp_s
        K4 = pcm_prop.k_s * pcm_prop.A_cond / (pcm_prop.z/2)

    elseif 0 <= L < 1
        Ln = L
        Mass_pcm = pcm_prop.ρ * pcm_prop.V
        cp = (pcm_prop.cp_l*(T[4]-pcm_prop.Tm1) + pcm_prop.Hlatent + pcm_prop.cp_s*((pcm_prop.Tm2-T[4])) ) / (pcm_prop.Tm2-pcm_prop.Tm1)
        # keff =   (pcm_prop.k_l*(T[4]-pcm_prop.Tm1) + pcm_prop.k_s*(pcm_prop.Tm2-T[4]) )/ (pcm_prop.Tm2-pcm_prop.Tm1)
        keff = pcm_prop.k_s + L * (pcm_prop.k_l-pcm_prop.k_s)
        K4 = keff * pcm_prop.A_cond / (pcm_prop.z/2)

    elseif L >= 1
        Ln = 1
        Mass_pcm = pcm_prop.ρ * pcm_prop.V
        cp = pcm_prop.cp_l
        K4 =  pcm_prop.k_l * pcm_prop.A_cond / (pcm_prop.z/2)

    end

    K4new = (1/Kpartial + 1/K4)^-1
    K[3,4] = K4new
    K[4,3] = K4new

    C[4] = Mass_pcm * cp

    H45new = (1/Hpartial + 1/K4)^-1
    H[4,5] = H45new
    H[5,4] = H45new

    # differential equations
    for i in 1:length(T)-1
        dT[i] = 1 ./ C[i] .* ( dot( H[i,:], (T .- T[i])) + dot( K[i,:], (T .- T[i])) + dot( R[i,:], (T.^4 .- T[i]^4) )  .+ Q[i] )
    end
    dT[6] = 0

end





end # module
