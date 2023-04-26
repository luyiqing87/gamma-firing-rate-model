# gamma-firing-rate-model
xpp code: XPP code for firing rate models.

    IN_rsv.ode: three-variable r-s-v firing rate for I-I network
    IN_rus.ode: three-variable r-u-s firing rate for I-I network
    IN_rs_delay.ode: two-variable r-s delay firing rate for I-I network
    EI_rsv.ode: six-variable r-s-v firing rate for E-I network
  
  
matlab code: MATLAB code for spiking network models and analysis of firing rate models.

    spknet_WB_II: Wang & Buzsaki I-I spiking network model
        WB_II.m: function that simulates the spiking network
        plot_rsv_II.m: script that calls WB_II to generate the time courses for r, s and v in Fig 2A
        
    spknet_HHlike_EI: an E-I spiking network model with Hodgkin-Huxley like single neuron 
        EI_network.m: function that simulates the spiking network
        ConEI.m: function that generates the connectivity matrix
        plot_rsv_EI.m: script that calls EI_network to generate the time courses for PING (set I_I=-1) and ING (set I_I=1) in Fig 1
        
    rate_rsv_II: three-variable r-s-v firing rate for I-I network
        rsv_parameters.mat: all parameters values for the r-s-v model as in Table 1
        rsv_solver.m: script that generates the time courses with MATLAB ODE solver ode45
        HB_I_finder.m: solve for the Hopf bifurcations with respect to the input drive I_I
        
    rate_rsdelay_II: two-variable r-s delay firing rate for I-I network
        parameters_rsdelay.mat: all parameters values for the r-s delay model as in Table 2
        rsdelay_solver.m: script that generates the time courses with MATLAB DDE solver dde23
        
    rate_rus_II: three-variable r-u-s firing rate for I-I network
        rus_parameters.mat: all parameters values for the r-u-s model as in Table 3
        rus_solver.m: script that generates the time courses with MATLAB ODE solver ode45
        
    rate_rsv_EI: six-variable r-s-v firing rate for E-I network
        rsv_EI_PINGpara.mat: all parameters values for the r-s-v model as in Table 4 for PING 
        rsv_EI_INGpara.mat: all parameters values for the r-s-v model as in Table 4 for ING except that I_I=0.5
        rsv_solver.m: script that generates the time courses with MATLAB ODE solver ode45
        nonlinsolver_6var_re_ri_plane.m: script that generates the rE, rI - SSRs in the rE-rI plane
