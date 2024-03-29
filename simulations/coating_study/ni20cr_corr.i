[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = 'EBSD_IC.txt'
  []
  uniform_refine = 2
  parallel_type = DISTRIBUTED
[]
[MultiApps]
  [grain_simulation]
    type = FullSolveMultiApp
    execute_on = INITIAL
    # positions = '0 0 0'
    input_files = 'ebsd_reader.i'
    # cli_args = 'GlobalParams/coating_thickness=${GlobalParams/coating_thickness}'
  []
[]

[Transfers]
  [to_corr_sim]
    source_variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1'
    variable =        'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1'
    type = MultiAppMeshFunctionTransfer
    from_multi_app = grain_simulation
  []
[]
[GlobalParams]
  op_num = 8
  var_name_base = gr
  derivative_order = 2
  enable_jit = false
  enable_ad_cache = false
  poly_grains = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11 gr12 gr13 gr14 gr15 gr16 gr17 gr18 gr19 gr20 gr21 gr22 gr23'
  int_width = 0.5

  K = ${fparse 8.6173324e-5 * 973 }
  x_Cr_metal = 0.2201
  x_Va = 2.2540e-07
  x_Ni_metal = ${fparse 1 - x_Va - x_Cr_metal }
  mu_Ni_alloy = ${fparse K*log(x_Ni_metal) - K*log(1.0 - 1.0*x_Ni_metal - 1.0*x_Cr_metal) - 1.729 }
  mu_Cr_alloy = ${fparse K*log(x_Cr_metal) - K*log(1.0 - 1.0*x_Ni_metal - 1.0*x_Cr_metal) - 1.717 }

  x_Cr_melt  = 25.0e-6
  x_Ni_melt  = 3.7872e-09 #3.758e-11
  E_F = 3.3607
  mu_Ni_salt = ${fparse K*log(x_Ni_melt) -5.5462 + 2*E_F }
  mu_Cr_salt = ${fparse K*log(x_Cr_melt) -6.8029 + 2*E_F }

  coating_thickness = %coating_thickness%
  x_Cr_coating = 1e-5
  x_Ni_coating = ${fparse 1 - x_Cr_coating - x_Va }
  mu_Ni_coating = ${fparse K*log(x_Ni_coating) - K*log(1.0 - x_Ni_coating - x_Cr_coating) - 1.729 }
  mu_Cr_coating = ${fparse K*log(x_Cr_coating) - K*log(1.0 - x_Ni_coating - x_Cr_coating) - 1.717 }
  # G_f_alloy = 5
  # G_f_coating = 7
[]

[Variables]
  [w_Ni]
  []
  [w_Cr]
  []
  [phi]
  []

  [eta1]
  []

  [./PolycrystalVariables]
  [../]
[]


[ICs]
  # # #Ni INITIAL CONDITIONS
  [w_Ni_global_initial]
    type = FunctionIC
    variable = 'w_Ni'
    function = 'mu_Ni_IC'
  []
  [w_Cr_global_initial]
    type = FunctionIC
    variable = 'w_Cr'
    function = 'mu_Cr_IC'
  []
[]

[Functions]
  [mu_Ni_IC]
    type = ParsedFunction
    vars = 'l L'
    vals = '${GlobalParams/int_width} ${GlobalParams/coating_thickness}'
    value = 'h1:=0.5*(1-tanh(2*(x-150.0)/l) );h2:=0.5*(1-tanh(2*(x-150-L)/l) );h1*${GlobalParams/mu_Ni_alloy} + (h2-h1)*${GlobalParams/mu_Ni_coating} + (1-h2)*${GlobalParams/mu_Ni_salt}'
  []
  [mu_Cr_IC]
    type = ParsedFunction
    vars = 'l L'
    vals = '${GlobalParams/int_width} ${GlobalParams/coating_thickness}'
    value = 'h1:=0.5*(1-tanh(2*(x-150.0)/0.5) );h2:=0.5*(1-tanh(2*(x-150-L)/l) );h1*${GlobalParams/mu_Cr_alloy} + (h2-h1)*${GlobalParams/mu_Cr_coating} + (1-h2)*${GlobalParams/mu_Cr_salt}'
  []
[]
[BCs]
  [phi_left]
    type = DirichletBC
    variable = phi
    value = 0
    boundary = 'left'
  []
  [w_Ni_right]
    type = DirichletBC
    variable = w_Ni
    value = '${GlobalParams/mu_Ni_salt}' #-0.45065
    boundary = 'right'
  []
  [w_Cr_right]
    type = FunctionDirichletBC
    variable = w_Cr
    function = '${GlobalParams/mu_Cr_salt}' #'-0.9699'
    boundary = 'right'
  []
[]

[AuxVariables]
  [bnds]
  []
  [E_x]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'INITIAL LINEAR TIMESTEP_END'
  []
  [electric_field]
    type = VariableGradientComponent
    variable = 'E_x'
    gradient_variable = 'phi'
    component = 'x'
  []
[]
[Kernels]
  #calculating phi
  [elec_chem_Ni]
    type = MatDiffusion
    variable = phi
    v = w_Ni
    diffusivity = 'chem_flux_Ni'
  []
  [elec_chem_Cr]
    type = MatDiffusion
    variable = phi
    v = w_Cr
    diffusivity = 'chem_flux_Cr'
  []
  [elec_phi]
    type = MatDiffusion
    variable = phi
    diffusivity = 'elec_flux'
  []

  # ##electro_diffusion
  [elec_diffusion_Ni]
    type = MatDiffusion
    variable = w_Ni
    v = phi
    diffusivity = elec_M_Ni
  []
  [elec_diffusion_Cr]
    type = MatDiffusion
    variable = w_Cr
    v = phi
    diffusivity = elec_M_Cr
  []
[]

[Modules]
  [PhaseField]
    [GrandPotential]
      switching_function_names = 'h_metal h_melt' #gr is metal, eta1 is electrolyte
      anisotropic = 'false false'

      chemical_potentials = 'w_Ni w_Cr'
      mobilities = 'M_Ni M_Cr'
      susceptibilities = 'chi_Ni chi_Cr'
      free_energies_w = 'c_Ni_metal c_Ni_melt c_Cr_metal c_Cr_melt'

      gamma_gr = gamma
      mobility_name_gr = L
      kappa_gr = kappa_gr
      free_energies_gr = 'omega_metal omega_melt'

      additional_ops = 'eta1'
      gamma_grxop = gamma
      mobility_name_op = L
      kappa_op = kappa
      free_energies_op = 'omega_metal omega_melt'

    []
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'gamma    gr_energy_sigma  interface_energy_sigma    interface_thickness_l   Va          pi           del_int   Na            xc  GB_width b l'
    prop_values = '1.5     6.803e6          1.429e7                    0.5                     1.109e-11  3.141592653  0.1       6.02214076e23 0.0 5e-04 0.8 0.01'
    outputs = 'nemesis'
  []
  [energy_constants]
    type = GenericConstantMaterial
    prop_names = 'kB            R      T     n  F          e        k_metal  k_melt   HF_H2  E0_F   E_F' #E0_F is electrode potential in Baes, E_F is pontential in the experimental salt
    prop_values = '8.6173324e-5 8.314  973  2  96485.33212 1.6e-19  1.0     1.0      1e-9   2.871  3.3607'
    outputs = 'nemesis'
  []

  [L]
    type = ParsedMaterial
    f_name = 'L'
    material_property_names = 'n F kB T Va interface_thickness_l pi del_int'
    constant_names = 'i0   m          D    J_to_eV    m3_to_um3 R'
    constant_expressions = '5.7 58.693e-3  8900 6.242e+18  1e18      8.314' #SI units
    function = '((del_int*i0*(m/D)^2)/(pi*n*F*R*T*interface_thickness_l*1e-6) )*(m3_to_um3/J_to_eV)'
  []

  [E0_Ni_metal]
    type = ParsedMaterial
    material_property_names = 'T F'
    function = '(-5179.159 + 117.854*T - 22.096*T*log(T) - (4.8407e-3)*T^2)/F'
    f_name = 'E0_Ni_metal'

  []
  [G_xs]
    type = ParsedMaterial
    f_name = 'Gxs'
    material_property_names = 'T F'
    constant_names = 'H_xs S_xs'
    constant_expressions = '-1.56448695e+04 -1.56011217'
    function = 'H_xs - S_xs*T'
  []
  [E0_Va_metal]
    type = ParsedMaterial
    material_property_names = 'T kB'
    constant_names = 'H0_f S0_f'
    constant_expressions = '1.56 3.3'
    function = 'H0_f - S0_f*kB*T'
    f_name = 'E0_Va_metal'
  []
  [E0_Cr_metal]
    type = ParsedMaterial
    material_property_names = 'T F Gxs'
    function = '(-1572.94 + 157.643*T - 26.908*T*log(T) + 1.89435e-3*T^2 - 1.47721e-6*T^3 + 139250/T '
               '+ Gxs)/F'
    f_name = 'E0_Cr_metal'
  []

  [E0_Ni_melt]
    type = ParsedMaterial
    f_name = 'E0_Ni_melt'
    constant_names = 'E0_NiF2'
    constant_expressions = '0.2973' #'0.473'
    material_property_names = 'E0_Ni_metal kB T n HF_H2 E0_F E_F'
    function = 'E0_Ni_metal + n*(E0_NiF2 + (E_F-E0_F))'
  []
  [E0_Cr_Cr]
    type = ParsedMaterial
    f_name = 'E0_Cr_Cr'
    material_property_names = 'T F'
    function = '(-8856.94 + 157.48*T - 26.908*T*log(T) + 1.89435e-3*T^2 - 1.47721e-6*T^3 + '
               '139250/T)/F'
  []
  [E0_Cr_melt]
    type = ParsedMaterial
    f_name = 'E0_Cr_melt'
    constant_names = 'E0_CrF2'
    constant_expressions = '-0.39'
    material_property_names = 'E0_Cr_Cr kB T n HF_H2 E0_F E_F'
    function = 'E0_Cr_Cr + n*(E0_CrF2 +  (E_F-E0_F))'
  []
  [E0_Va_melt]
    type = GenericConstantMaterial
    prop_names = 'E0_Va_melt'
    prop_values = '0.0'
  []
  #PARAMETERS
  [kappa] #assume that three interfaces having the same interfacial energy and thickness
    type = ParsedMaterial
    f_name = kappa
    material_property_names = 'interface_energy_sigma interface_thickness_l'
    function = '3*interface_energy_sigma*interface_thickness_l/4'
  []
  [kappa_gr] #assume that three interfaces having the same interfacial energy and thickness
    type = ParsedMaterial
    f_name = kappa_gr
    material_property_names = 'gr_energy_sigma interface_thickness_l'
    function = '3*gr_energy_sigma*interface_thickness_l/4'
  []
  [m]
    type = ParsedMaterial
    f_name = mu
    material_property_names = 'interface_energy_sigma interface_thickness_l'
    function = '6*interface_energy_sigma/interface_thickness_l'
  []
  #SWITCHING FUNCTIONS
  [switch_melt]
    type = DerivativeParsedMaterial
    f_name = 'h_melt'
    args = 'eta1 ${GlobalParams/poly_grains}'
    material_property_names = 'pi del_int'
    function = '0.5*(1 + tanh( 2*pi*(eta1-0.5)/del_int ) )'
    output_properties = 'h_melt'
  []
  [switch_metal]
    type = DerivativeParsedMaterial
    f_name = 'h_metal'
    args = 'eta1 ${GlobalParams/poly_grains}'
    material_property_names = 'pi del_int'
    function = '0.5*(1 - tanh( 2*pi*(eta1 -0.5)/del_int ) )'
    output_properties = 'h_metal'
  []
  [h_gb]
    type = ParsedMaterial
    args = '${GlobalParams/poly_grains} bnds'
    f_name = 'h_gb'
    material_property_names = 'pi del_int b l'
    function = '0.5-0.5*tanh(pi*(bnds-b)/l)'
    outputs = nemesis
    output_properties = 'h_gb'
  []
  [c_Ni_metal]
    type = DerivativeParsedMaterial
    f_name = "c_Ni_metal"
    args = 'w_Ni w_Cr'
    material_property_names = 'kB T k_metal E0_Ni_metal E0_Va_metal  E0_Cr_metal'
    function = 'exp( (w_Ni - (E0_Ni_metal - E0_Va_metal))/kB/T/k_metal )/( 1 + exp( (w_Ni - '
               '(E0_Ni_metal - E0_Va_metal))/kB/T/k_metal ) + exp( (w_Cr - (E0_Cr_metal - '
               'E0_Va_metal))/kB/T/k_metal ))'
    output_properties = 'c_Ni_metal'
  []
  [c_Ni_melt]
    type = DerivativeParsedMaterial
    f_name = "c_Ni_melt"
    args = 'w_Ni w_Cr'
    material_property_names = 'kB T k_melt E0_Ni_melt E0_Va_melt E0_Cr_melt'
    function = 'exp( (w_Ni - (E0_Ni_melt) - k_melt*kB*T)/kB/T/k_melt )'
    output_properties = 'c_Ni_melt'
  []
  [c_Cr_metal]
    type = DerivativeParsedMaterial
    f_name = "c_Cr_metal"
    args = 'w_Ni w_Cr'
    material_property_names = 'kB T k_metal E0_Ni_metal E0_Cr_metal E0_Va_metal'
    function = 'exp( (w_Cr - (E0_Cr_metal - E0_Va_metal))/kB/T/k_metal )/( 1 + exp( (w_Ni - '
               '(E0_Ni_metal - E0_Va_metal))/kB/T/k_metal ) + exp( (w_Cr - (E0_Cr_metal - '
               'E0_Va_metal))/kB/T/k_metal ))'
    output_properties = 'c_Cr_metal'

  []
  [c_Cr_melt]
    type = DerivativeParsedMaterial
    f_name = "c_Cr_melt"
    args = 'w_Ni w_Cr'
    material_property_names = 'kB T k_melt E0_Ni_melt E0_Va_melt E0_Cr_melt'
    function = 'exp( (w_Cr - (E0_Cr_melt - E0_Va_melt) - k_melt*kB*T)/kB/T/k_melt )'
    #outputs = exodus
    output_properties = 'c_Cr_melt'

  []
  #Calculate global concentrations using switching FUNCTIONS
  [c_Ni]
    type = DerivativeParsedMaterial
    f_name = 'c_Ni'
    material_property_names = 'c_Ni_metal c_Ni_melt h_metal h_melt'
    function = 'c_Ni_metal*h_metal + c_Ni_melt*(1-h_metal)'
    output_properties = 'c_Ni'
    outputs = nemesis
  []
  [c_Cr]
    type = DerivativeParsedMaterial
    f_name = 'c_Cr'
    material_property_names = 'c_Cr_metal c_Cr_melt h_metal h_melt'
    function = 'c_Cr_metal*h_metal + c_Cr_melt*(1-h_metal)'
    output_properties = 'c_Cr'
    outputs = nemesis
  []
  [c_Va]
    type = DerivativeParsedMaterial
    f_name = 'c_Va'
    material_property_names = 'c_Ni c_Cr'
    function = '1-c_Ni-c_Cr'
    outputs = nemesis
    output_properties = 'c_Va'
  []
  [del_c_Ni]
    type = DerivativeParsedMaterial
    f_name = 'del_c_Ni'
    args = 'w_Ni w_Cr ${GlobalParams/poly_grains} eta1 '
    material_property_names = 'kB T k_metal E0_Ni_metal E0_Va_metal  E0_Cr_metal k_melt E0_Ni_melt '
                              'E0_Va_melt E0_Cr_melt'
    function = 'exp( (w_Ni - (E0_Ni_metal - E0_Va_metal))/kB/T/k_metal )/( 1 + exp( (w_Ni - '
               '(E0_Ni_metal - E0_Va_metal))/kB/T/k_metal ) + exp( (w_Cr - (E0_Cr_metal - '
               'E0_Va_metal))/kB/T/k_metal ))
     - exp( (w_Ni - (E0_Ni_melt-E0_Va_melt) - '
               'k_melt*kB*T )/kB/T/k_melt )'
  []
  [del_c_Cr]
    type = DerivativeParsedMaterial
    f_name = 'del_c_Cr'
    args = 'w_Ni w_Cr ${GlobalParams/poly_grains} eta1 '
    material_property_names = 'kB T k_metal E0_Ni_metal E0_Cr_metal E0_Va_metal k_melt E0_Ni_melt '
                              'E0_Va_melt E0_Cr_melt'
    function = 'exp( (w_Cr - (E0_Cr_metal - E0_Va_metal))/kB/T/k_metal )/( 1 + exp( (w_Ni - '
               '(E0_Ni_metal - E0_Va_metal))/kB/T/k_metal ) + exp( (w_Cr - (E0_Cr_metal - '
               'E0_Va_metal))/kB/T/k_metal ))
     - exp( (w_Cr - (E0_Cr_melt-E0_Va_melt)  - '
               'k_melt*kB*T)/kB/T/k_melt )'
  []
  #Derivative terms of free energy
  [omega_metal]
    type = DerivativeParsedMaterial
    f_name = 'omega_metal'
    args = '${GlobalParams/poly_grains} eta1 w_Ni w_Cr'
    material_property_names = 'kB T k_metal E0_Ni_metal E0_Va_metal  E0_Cr_metal Va'
    function = '(E0_Va_metal - kB*T*k_metal*log( 1 + exp( (w_Ni - (E0_Ni_metal - '
               'E0_Va_metal))/kB/T/k_metal ) + exp( (w_Cr - (E0_Cr_metal - '
               'E0_Va_metal))/kB/T/k_metal ) ) )/Va'
    derivative_order = 2
    outputs = nemesis
    output_properties = 'omega_metal'
  []
  [omega_melt]
    type = DerivativeParsedMaterial
    f_name = 'omega_melt'
    args = '${GlobalParams/poly_grains}  eta1 w_Ni w_Cr'
    material_property_names = 'kB T k_melt E0_Ni_melt E0_Va_melt E0_Cr_melt Va'
    function = '(E0_Va_melt - kB*T*k_melt*(exp( (w_Ni - (E0_Ni_melt - E0_Va_melt) - '
               'k_melt*kB*T)/kB/T/k_melt ) + exp( (w_Cr - (E0_Cr_melt - E0_Va_melt) - '
               'k_melt*kB*T)/kB/T/k_melt ) ) )/Va'
    outputs = nemesis
    output_properties = 'omega_melt'
    derivative_order = 2
  []
  [omega_chem]
    type = DerivativeParsedMaterial
    f_name = omega_chem
    args = '${GlobalParams/poly_grains} eta1 w_Ni w_Cr'
    material_property_names = 'h_metal omega_metal omega_melt Va'
    function = '(omega_metal-omega_melt)*Va'
    outputs = nemesis
    output_properties = 'omega_chem'
  []
  [susceptibility_Ni]
    type = DerivativeParsedMaterial
    f_name = 'chi_Ni'
    args = '${GlobalParams/poly_grains} eta1  w_Ni w_Cr'
    material_property_names = 'c_Ni_metal(w_Ni,w_Cr) c_Ni_melt(w_Ni,w_Cr) h_metal '
                              'chi_Ni_metal:=D[c_Ni_metal,w_Ni] chi_Ni_melt:=D[c_Ni_melt,w_Ni]'
    function = 'chi_Ni_metal*h_metal + chi_Ni_melt*(1-h_metal)'
    outputs = nemesis
    output_properties = 'chi_Ni'
  []
  [D_Ni_V]
    type = ParsedMaterial
    f_name = 'D_Ni_V'
    material_property_names = 'T R c_Va'
    constant_names = 'D0_Ni_V       E0_Ni_V     cal_to_J' ##https://aip.scitation.org/doi/pdf/10.1063/1.1703047
    constant_expressions = '1.9e8         66800       4.184'
    function = 'D0_Ni_V*exp(-E0_Ni_V*cal_to_J/R/T)*(c_Va/2.254e-7)'
    outputs = nemesis
    output_properties = 'D_Ni_V'
  []
  [D_Ni_GB]
    type = ParsedMaterial
    f_name = 'D_Ni_GB'
    material_property_names = 'T R GB_width c_Va'
    constant_names = 'D0_Ni_GB      E0_Ni_GB     cal_to_J' ##https://aip.scitation.org/doi/pdf/10.1063/1.1703047
    constant_expressions = '0.07e8        27400        4.184'
    function = 'D0_Ni_GB*exp(-E0_Ni_GB*cal_to_J/R/T)*(c_Va/2.254e-7)'
    outputs = nemesis
    output_properties = 'D_Ni_GB'
  []
  [D_Ni]
    type = ParsedMaterial
    f_name = 'D_Ni'
    material_property_names = 'D_Ni_V D_Ni_GB h_metal h_gb  interface_thickness_l GB_width'
    function = '(D_Ni_V +   1.0*(D_Ni_GB-D_Ni_V )*(GB_width/interface_thickness_l)*h_gb)*h_metal' # + 500*h_melt'
    outputs = nemesis
    output_properties = 'D_Ni'
  []
  [mobility_Ni]
    type = DerivativeParsedMaterial
    f_name = M_Ni
    args = '${GlobalParams/poly_grains} eta1 w_Ni w_Cr'
    material_property_names = 'chi_Ni h_metal h_melt h_gb D_Ni_V D_Ni_GB interface_thickness_l '
                              'GB_width'
    function = '((D_Ni_V +   1.0*(D_Ni_GB-D_Ni_V )*(GB_width/interface_thickness_l)*h_gb)*h_metal + '
               '500*h_melt )*chi_Ni'
    outputs = nemesis
    output_properties = 'M_Ni'
    derivative_order = 2
  []

  [susceptibility_Cr]
    type = DerivativeParsedMaterial
    f_name = 'chi_Cr'
    args = '${GlobalParams/poly_grains} eta1 w_Ni w_Cr'
    material_property_names = 'c_Cr_metal(w_Ni,w_Cr) c_Cr_melt(w_Ni,w_Cr) h_metal '
                              'chi_Cr_metal:=D[c_Cr_metal,w_Cr] chi_Cr_melt:=D[c_Cr_melt,w_Cr]'
    function = 'chi_Cr_metal*h_metal + chi_Cr_melt*(1-h_metal)'
    outputs = nemesis
    output_properties = 'chi_Cr'
    #outputs = exodus
  []
  [D_Cr_V]
    type = ParsedMaterial
    f_name = 'D_Cr_V'
    material_property_names = 'T R xc c_Va'
    function = '(exp(-(84702800778885824000*xc + '
               '306299812110263875)/(8927089524736*T))*exp((5885692387260437*xc)/1099511627776)*exp(6'
               '99261826406127/35184372088832))*(c_Va/2.254e-7)'
    outputs = nemesis
    output_properties = 'D_Cr_V'
  []
  [D_Cr_GB]
    type = ParsedMaterial
    f_name = 'D_Cr_GB'
    material_property_names = 'T R xc interface_thickness_l GB_width c_Va'
    function = '(exp((209787039506469*xc)/17179869184)*exp(2393065153853383/140737488355328)*exp(-(85'
               '3955374395308288000*xc + '
               '841876566076058375)/(35708358098944*T)))*(c_Va/2.254e-7)/GB_width'
    outputs = nemesis
    output_properties = 'D_Cr_GB'
  []
  [D_Cr]
    type = ParsedMaterial
    f_name = 'D_Cr'
    material_property_names = 'D_Cr_V D_Cr_GB h_metal h_gb interface_thickness_l GB_width'
    function = '(D_Cr_V +   1.0*(D_Cr_GB-D_Cr_V )*(GB_width/interface_thickness_l)*h_gb)*h_metal'
    outputs = nemesis
    output_properties = 'D_Cr'
  []
  [mobility_Cr]
    type = DerivativeParsedMaterial
    f_name = M_Cr
    args = '${GlobalParams/poly_grains} eta1  w_Ni w_Cr'
    material_property_names = 'chi_Cr h_metal h_melt h_gb interface_thickness_l GB_width D_Cr_V '
                              'D_Cr_GB'
    function = '((D_Cr_V +   1.0*( D_Cr_GB - D_Cr_V )*(GB_width/interface_thickness_l)*h_gb)*h_metal '
               '+ 500*h_melt )*chi_Cr'
    outputs = nemesis
    output_properties = 'M_Cr'
    derivative_order = 2
  []
  [M_H]
    type = DerivativeParsedMaterial
    args = '${GlobalParams/poly_grains} eta1  w_Ni w_Cr'
    f_name = 'M_H'
    material_property_names = 'h_metal h_melt kB T k_melt R E0_F E_F'

    function = 'c_H:= exp(-(E_F-E0_F)/kB/T);(1e-7*h_metal + 1000*h_melt )*c_H/kB/T/k_melt'
    derivative_order = 2
  []
  ##For Nichrome,rho = 1.5e-6 ohm-m (T=293 K), alpha = 0.4e-3. At T = 973 K, rho_T = 1.5e-6*(1+0.4e-3*(973-293) ) = 1.90800e-6 ohm-m
  ##Then sigma_e = 1/rho_T = = 5.24e5 S/m = 5.24e-1 S/micrometer
    [elec_conductivity]
      type = ParsedMaterial
      f_name = 'sigma_e'
      function = '0.524' #S/um
      outputs = 'nemesis'
      output_properties = 'sigma_e'
    []

  [M_e]
    type = DerivativeParsedMaterial
    args = 'w_Ni ${GlobalParams/poly_grains} eta1 '
    f_name = 'M_e'
    material_property_names = 'h_metal h_melt kB T k_metal Va sigma_e e'
    function = '( (sigma_e*Va/2/e)*h_metal + 1e-12*h_melt )/kB/T/k_metal'
    derivative_order = 2
  []
  #Computing phi
  [flux_Ni]
    type = DerivativeParsedMaterial
    args = '${GlobalParams/poly_grains} eta1  w_Ni phi'
    f_name = 'chem_flux_Ni'
    material_property_names = '   M_Ni M_H h_metal h_melt'
    function = 'h_melt*M_Ni*2'
    derivative_order = 2
  []
  [flux_Cr]
    type = DerivativeParsedMaterial
    args = '${GlobalParams/poly_grains} eta1 w_Cr phi'
    f_name = 'chem_flux_Cr'
    material_property_names = '   M_Cr M_H h_metal h_melt'
    function = 'h_melt*M_Cr*2'
    output_properties = 'chem_flux_Cr'
    derivative_order = 2
  []
  [grad_phi_coeff]
    type = DerivativeParsedMaterial
    f_name = 'elec_flux'
    args = '${GlobalParams/poly_grains} eta1  w_Ni phi'
    material_property_names = '   M_Ni M_Cr M_H M_e h_metal h_melt'
    function = 'h_melt*(4*M_Ni + 4*M_Cr + M_H) + (h_metal)*M_e'
    output_properties = 'elec_flux'
    derivative_order = 2
  []

  #electro_diffusion
  [elec_M_Ni]
    type = DerivativeParsedMaterial
    f_name = 'elec_M_Ni'
    args = 'w_Ni ${GlobalParams/poly_grains} eta1  phi'
    material_property_names = 'h_metal h_melt M_Ni'
    function = 'h_melt*M_Ni*2'
    derivative_order = 2
  []
  [elec_M_Cr]
    type = DerivativeParsedMaterial
    f_name = 'elec_M_Cr'
    args = 'w_Cr ${GlobalParams/poly_grains} eta1  phi'
    material_property_names = 'h_metal h_melt M_Cr'
    function = 'h_melt*M_Cr*2'
    derivative_order = 2
  []

  #Postprocessor materials
  [metal_thickness]
    type = ParsedMaterial
    f_name = 'metal_thickness'
    material_property_names = 'h_metal'
    function = 'if(h_metal>0.5,1,0)'
  []
  [Cr_mass]
    type = ParsedMaterial
    f_name = 'Cr_mass'
    material_property_names = 'c_Cr h_metal metal_thickness Va'
    function = '(h_metal*c_Cr/Va)*51.99*1e8*1000/6.022e23' #g/um3 -> 1e12 g/cm3, Mass loss/Area => g/um2 -> 1e8 g/cm2, *1000 -> mg/cm2
  []
  [Ni_mass]
    type = ParsedMaterial
    f_name = 'Ni_mass'
    material_property_names = 'c_Ni h_metal metal_thickness Va'
    function = '(h_metal*c_Ni/Va)*58.69*1e8*1000/6.022e23' #g/um3 -> 1e12 g/cm3, Mass loss/Area => g/um2 -> 1e8 g/cm2, *1000 -> mg/cm2
  []
  [total_mass]
    type = ParsedMaterial
    f_name = 'total_mass'
    material_property_names = 'c_Cr_metal c_Ni_metal h_metal metal_thickness Va'
    function = '(h_metal/Va)*(c_Cr_metal*51.99 +c_Ni_metal*58.69)*1e8*1000/6.022e23' #g/um3 -> 1e12 g/cm3, Mass loss/Area => g/um2 -> 1e8 g/cm2, *1000 -> mg/cm2
  []

[]
[Postprocessors]
  [elapsed]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
    outputs = 'csv'
  []
  [Cr_metal_mass_total]
    type = ElementIntegralMaterialProperty
    mat_prop = 'Cr_mass'
    outputs = 'csv'
  []
  [Ni_metal_mass_total]
    type = ElementIntegralMaterialProperty
    mat_prop = 'Ni_mass'
    outputs = 'csv'
  []
  [total_alloy_mass]
    type = ElementIntegralMaterialProperty
    mat_prop = 'total_mass'
    outputs = 'csv'
  []
  [total_Ni]
    type = ElementIntegralMaterialProperty
    mat_prop = 'c_Ni'
    outputs = 'csv'
  []
  [total_Cr]
    type = ElementIntegralMaterialProperty #ElementIntegralMaterialProperty
    mat_prop = 'c_Cr'
    outputs = 'csv'
  []
  [Num_DOFs]
    type = NumDOFs
    execute_on = 'initial timestep_end'
    outputs = 'csv'
    system = 'ALL'
  []
  [Avg_D_Cr]
    type = ElementAverageMaterialProperty
    mat_prop = 'D_Cr'
    outputs = 'csv'
  []
  [total_metal_phase]
    type = ElementIntegralMaterialProperty
    mat_prop = 'h_metal'
  []
  [corrosion_current]
    type = SideFluxAverage
    variable = phi
    diffusivity = 'sigma_e'
    boundary = left
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK #
  scheme = bdf2
  # petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  # petsc_options_value = 'asm      31                  preonly       lu           2'
  # petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold'
  # petsc_options_value = 'hypre    boomeramg      0.7'
  petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl  -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor'
  petsc_options_value = 'hypre    boomeramg      0.7                                   4                          5  25  HMIS  ext+i  2  0.3'
  l_max_its = 30
  nl_max_its = 10
  l_tol = 1e-3
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
  dtmin = 1e-7
  end_time = 3.6e6
  automatic_scaling = true
  dtmax = 1e4
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    growth_factor = 1.25
    cutback_factor = 0.8
    iteration_window = 2
    optimal_iterations = 7
  []
  # num_steps = 0
[]

[Outputs]
  [nemesis]
    type = Nemesis
    execute_on = 'initial timestep_end final'
    interval = 10
  []
  csv = true
  perf_graph = true
  file_base = 'ni20_corr/ni20_corr'
[]
