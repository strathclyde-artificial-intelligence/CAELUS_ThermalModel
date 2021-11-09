#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 14:38:17 2021

@author: gianlucafilippi
paper: https://www.readcube.com/articles/10.3389/fmech.2019.00029
"""


def model_pcm(temperature, params, C, K, H):
    # PCM interpolation
    K4 = []
    norm_temperature = (temperature - params["PCM_t1"]) / (params["PCM_t2"] - params["PCM_t1"])

    if temperature < params["PCM_t1"]:
        cp = params["PCM_cps"]

        Ln = 0
        Mass_pcm = params["PCM_rho"] * params["PCM_V"]
        K4 = params["PCM_ks"] * params["PCM_A_cond"] / (params["PCM_z"] / 2)

    elif params["PCM_t1"] <= temperature < params["PCM_t2"]:
        Ln = norm_temperature
        Mass_pcm = params["PCM_rho"] * params["PCM_V"]

        cp = (params["PCM_cps"] * (params["PCM_t2"] - temperature)
              + params["PCM_cpl"] * (temperature - params["PCM_t1"])
              + params["PCM_h"]) \
            / (params["PCM_t2"] - params["PCM_t1"])

        keff = params["PCM_ks"] + norm_temperature * (params["PCM_kl"] - params["PCM_ks"])
        K4 = keff * params["PCM_A_cond"] / (params["PCM_z"] / 2)

    else:
        Ln = 1
        cp = params["PCM_cpl"]
        Mass_pcm = params["PCM_rho"] * params["PCM_V"]
        K4 = params["PCM_kl"] * params["PCM_A_cond"] / (params["PCM_z"] / 2)

    Kpartial = K[3, 2]
    Hpartial = H[3, 4]

    # update 
    K4new = pow((1 / Kpartial + 1 / K4), -1)
    K[2, 3] = K4new
    K[3, 2] = K4new

    C[3] = Mass_pcm * cp

    H45new = pow((1 / Hpartial + 1 / K4), -1)
    H[3, 4] = H45new
    H[4, 3] = H45new

    return K, H
