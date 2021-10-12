import numpy as np

def thermal_model(t, T, PCM_model, atmos_model, params):
    """
    Explain me!
    """ 
    dT = [0, 0, 0, 0, 0, 0]
    
    #self.params, M, self.atmos_model = params
    
    # C = ode_matrices.get_C()
    # K = ode_matrices.get_K()
    # H = ode_matrices.get_H()
    # R = ode_matrices.get_R()
    # Q = ode_matrices.get_Q()
    
    
    C = params["matrix_C"]
    K = params["matrix_K"]
    H = params["matrix_H"]
    R = params["matrix_R"]
    Q = params["matrix_Q"]


    # ambient temperature from Atmosferic model
    T[5] = atmos_model(t)
    
    
    K, H = PCM_model(T, params, C, K, H)
    
    
    # Kpartial = K[3,2]
    # Hpartial = H[3,4]

    # # C, K, H, R, Q, Kpartial, Hpartial, pcm_prop, atmmodel = params





    # # PCM interpolation
    # K4 = []
    # L = (T[3]-params["PCM_t1"]) / (params["PCM_t2"]-params["PCM_t1"])
    # # PCM properties
    # if L < 0:
    #     Ln = 0
    #     Mass_pcm = params["PCM_rho"] * params["PCM_V"]
    #     cp = params["PCM_cps"]
    #     K4 = params["PCM_ks"] * params["PCM_A_cond"] / (params["PCM_z"]/2)

    # elif (0 <= L < 1):
    #     Ln = L
    #     Mass_pcm = params["PCM_rho"] * params["PCM_V"]
    #     cp = (params["PCM_cpl"] * (T[3] - params["PCM_t1"]) + params["PCM_h"] + params["PCM_cps"] * ((params["PCM_t2"] - T[3])) ) / (params["PCM_t2"] - params["PCM_t1"])
    #     # keff =   (pcm_prop.k_l*(T[4]-pcm_prop.Tm1) + pcm_prop.k_s*(pcm_prop.Tm2-T[4]) )/ (pcm_prop.Tm2-pcm_prop.Tm1)
    #     keff = params["PCM_ks"] + L * (params["PCM_kl"] - params["PCM_ks"])
    #     K4 = keff * params["PCM_A_cond"] / (params["PCM_z"]/2)

    # else:
    #     Ln = 1
    #     Mass_pcm = params["PCM_rho"] * params["PCM_V"]
    #     cp = params["PCM_cpl"]
    #     K4 = params["PCM_kl"] * params["PCM_A_cond"] / (params["PCM_z"]/2)

    

    # K4new = pow((1/Kpartial + 1/K4), -1)
    # K[2,3] = K4new
    # K[3,2] = K4new

    # C[3] = Mass_pcm * cp

    # H45new = pow((1/Hpartial + 1/K4), -1)
    # H[3,4] = H45new
    # H[4,3] = H45new

    TT_row = np.array([ [T[0]]*6, [T[1]]*6, [T[2]]*6, [T[3]]*6, [T[4]]*6, [T[5]]*6 ])
    TT_col = TT_row.transpose()
    
    TT = TT_row - TT_col
    TTR = pow(TT_row, 4) - pow(TT_col, 4)

    dT = 1./C * (np.diag(np.dot(H, TT)) + np.diag(np.dot(K, TT)) + np.diag(np.dot(R, TTR)) + Q)
    dT[5] = 0

    return dT
