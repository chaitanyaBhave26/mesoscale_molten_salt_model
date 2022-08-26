[Mesh]
  [gen]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 160 ##2000 #40 #200
    xmin = -1e-5
    xmax = 159.53 ##999.99999
    ny = 300
    ymin = -1e-5
    ymax = 299.99999 ##24.99999 #100
    # elem_type = QUAD4
  []

  uniform_refine = 3#9 #6#8
  parallel_type = DISTRIBUTED
[]
# [Mesh]
#   type = FileMesh
#   file = 'multigrain_out.e'
#   uniform_refine = 3
# []
[GlobalParams]
  op_num = 8
  var_name_base = gr
  # outputs = nemesis
  derivative_order = 2
  enable_jit = false
  enable_ad_cache = false
  # int_width = 10#0.5
  # grain_num = 1#2
  # numbub=1
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

  [gr0]
  []
  [gr1]
  []
  [gr2]
  []
  [gr3]
  []
  [gr4]
  []
  [gr5]
  []
  [gr6]
  []
  [gr7]
  []
  # [gr8]
  # []
  # [gr9]
  # []
  # [gr10]
  # []
  # [gr11]
  # []
  # [gr12]
  # []
  # [gr13]
  # []
  # [gr14]
  # []

[]

[UserObjects]
  [eta_init]
    type = SolutionUserObject
    # mesh = 'ebsd_reader/ebsd_reader.e-s042' #'multigrain/multigrain.e-s087'
    mesh = 'ebsd_reader/ebsd_reader.e'#'multigrain/multigrain.e-s087'
    # mesh = 'ebsd_reader1/ebsd_reader.e'#'multigrain/multigrain.e-s087'
    # mesh = '/blue/michael.tonks/share/ebsd_reader.e-s012'#'multigrain/multigrain.e-s087'
    system_variables = 'eta1 gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11'# gr12 gr13 gr14'
    # mesh = '2_grain_ic.e-s016'
    # system_variables = 'eta_melt gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    timestep = 'LATEST'
    # execute_on = 'INITIAL'
  []
  # [eta_init]
  #   type = SolutionUserObject
  #   mesh = 'multigrain/multigrain.e-s168'#'multigrain/multigrain.e-s087'
  #   # system_variables = 'eta1 gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
  #   # mesh = '2_grain_ic.e-s016'
  #   system_variables = 'eta_melt gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
  #   timestep = 'LATEST'
  # []
[]
[Functions]
  [gr0_val]
    type = SolutionFunction
    solution = eta_init
    from_variable = 'gr0'
  []
  [gr1_val]
    type = SolutionFunction
    solution = eta_init
    from_variable = 'gr1'
  []
  [gr2_val]
    type = SolutionFunction
    solution = eta_init
    from_variable = 'gr2'
  []
  [gr3_val]
    type = SolutionFunction
    solution = eta_init
    from_variable = 'gr3'
  []
  [gr4_val]
    type = SolutionFunction
    solution = eta_init
    from_variable = 'gr4'
  []
  [gr5_val]
    type = SolutionFunction
    solution = eta_init
    from_variable = 'gr5'
  []
  [gr6_val]
    type = SolutionFunction
    solution = eta_init
    from_variable = 'gr6'
  []
  [gr7_val]
    type = SolutionFunction
    solution = eta_init
    from_variable = 'gr7'
  []
  # [gr8_val]
  #   type = SolutionFunction
  #   solution = eta_init
  #   from_variable = 'gr8'
  # []
  # [gr9_val]
  #   type = SolutionFunction
  #   solution = eta_init
  #   from_variable = 'gr9'
  # []
  # [gr10_val]
  #   type = SolutionFunction
  #   solution = eta_init
  #   from_variable = 'gr10'
  # []
  # [gr11_val]
  #   type = SolutionFunction
  #   solution = eta_init
  #   from_variable = 'gr11'
  # []
  # [gr12_val]
  #   type = SolutionFunction
  #   solution = eta_init
  #   from_variable = 'gr12'
  # []
  # [gr13_val]
  #   type = SolutionFunction
  #   solution = eta_init
  #   from_variable = 'gr13'
  # []
  # [gr14_val]
  #   type = SolutionFunction
  #   solution = eta_init
  #   from_variable = 'gr14'
  # []
  [eta1_val]
    type = SolutionFunction
    solution = eta_init
    from_variable = 'eta1'
  []
[]
[ICs]
  [gr0_init]
    type = FunctionIC
    function = 'gr0_val'
    variable = 'gr0'
  []
  [gr1_init]
    type = FunctionIC
    function = 'gr1_val'
    variable = gr1
  []
  [gr2_init]
    type = FunctionIC
    function = 'gr2_val'
    variable = gr2
  []
  [gr3_init]
    type = FunctionIC
    function = 'gr3_val'
    variable = gr3
  []
  [gr4_init]
    type = FunctionIC
    function = 'gr4_val'
    variable = gr4
  []
  [gr5_init]
    type = FunctionIC
    function = 'gr5_val'
    variable = gr5
  []
  [gr6_init]
    type = FunctionIC
    function = 'gr6_val'
    variable = gr6
  []
  [gr7_init]
    type = FunctionIC
    function = 'gr7_val'
    variable = gr7
  []
  # [gr8_init]
  #   type = FunctionIC
  #   function = 'gr8_val'
  #   variable = 'gr8'
  # []
  # [gr9_init]
  #   type = FunctionIC
  #   function = 'gr9_val'
  #   variable = gr9
  # []
  # [gr10_init]
  #   type = FunctionIC
  #   function = 'gr10_val'
  #   variable = gr10
  # []
  # [gr11_init]
  #   type = FunctionIC
  #   function = 'gr11_val'
  #   variable = gr11
  # []
  # [gr12_init]
  #   type = FunctionIC
  #   function = 'gr12_val'
  #   variable = gr12
  # []
  # [gr13_init]
  #   type = FunctionIC
  #   function = 'gr13_val'
  #   variable = gr13
  # []
  # [gr14_init]
  #   type = FunctionIC
  #   function = 'gr14_val'
  #   variable = gr14
  # []
  [eta1_init]
    type = FunctionIC
    function = 'eta1_val'
    variable = eta1
  []
  # #Ni INITIAL CONDITIONS
  [c_global_inital]
    type = BoundingBoxIC
    variable = 'w_Ni'
    x1 = -20
    y1 = -20
    x2 = 300
    y2 = 520
    inside = -0.45065#-0.4667
    outside = -0.4667 #1.936039
  []
  [c_Cr_global_inital]
    type = FunctionIC
    variable = 'w_Cr'
    # function = '0.5-tanh(2*(x-320)/5)'

    # function = '-0.56089'#'h1:=0.5-0.5*tanh(2*(x-320)/5);-0.56089*h1 + (1-h1)*(-0.9699)'
    # function = 'if(x<300,-0.56089,-0.5608 - (-0.5608 + 0.9699)*(x-300)/800 )'
    # function = 'if(x<300,-0.56089,-0.5608 - (-0.5608 + 0.9699)*(x-300)/200 )'
    # function = 'if(x<300,-0.56089,-0.5608 - (-0.5608 + 0.7699)*(x-300)/200 )'
    ##function = 'if(x<300,-0.56089,-0.5608 - (-0.5608 + 0.9699)*(x-300)/200 )'
    function = '-0.6755'
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
    variable = w_Cr
    value = -0.4667
    boundary = 'right'
  []
  [w_Cr_right]
    type = FunctionDirichletBC
    variable = w_Cr
    function = '-0.9699'#'if(t<100,1.9764 - (1.9764-0.83884)/sqrt(1.5e7/t),0.83884)'
    # function = 'if(t<100,1.9764 - (1.9764-0.83884)/sqrt(1.5e7/t),0.83884)'
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
  [eta2]
    order = FIRST
    family = LAGRANGE
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
  [./PhaseField]
    [./GrandPotential]
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

    [../]
  [../]
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'gamma    gr_energy_sigma  interface_energy_sigma    interface_thickness_l  Va pi '
                 'del_int Na xc GB_width b l'
    prop_values = '1.5     10                   5.1e-5                          0.50         1.1131e-11 '
                  ' 3.141592653 0.1 6.02214076e23 0.0 5e-04 0.8 0.01'
    outputs = 'nemesis'
  []
  [energy_constants]
    type = GenericConstantMaterial
    prop_names = 'kB            R      T     n  F             k_metal  k_melt   HF_H2  E0_F   E_F' #E0_F is electrode potential in Baes, E_F is pontential in the experimental salt
    prop_values = '8.6173324e-5 8.314  973  2  96485.33212   1.0     1.0      1e-9   2.871  3.3607' #12e-4
    outputs = 'nemesis'
  []

  [L]
    type = ParsedMaterial
    f_name = 'L'
    material_property_names = 'n F kB T Va interface_thickness_l pi del_int'
    constant_names = 'i0   m          D    J_to_eV    m3_to_um3 R'
    constant_expressions = '5.7 58.693e-3  8900 6.242e+18  1e18      8.314' #SI units
    function = '((del_int*i0*(m/D)^2)/(pi*n*F*R*T*interface_thickness_l*1e-6) '
               ')*(m3_to_um3/J_to_eV)/Va' #/Va to normalize L for grand potential in eV/atom rather than eV/um3
    #outputs = exodus
  []

  [E0_Ni_metal]
    type = ParsedMaterial
    material_property_names = 'T F'
    function = '(-5179.159 + 117.854*T - 22.096*T*log(T) - (4.8407e-3)*T^2)/F'
    f_name = 'E0_Ni_metal'
    # outputs = exodus
  []
  [G_xs]
    type = ParsedMaterial
    f_name = 'Gxs'
    material_property_names = 'T F'
    constant_names = 'H_xs S_xs'
    constant_expressions = '-1.56448695e+04 -1.56011217'
    function = 'H_xs - S_xs*T' #'H_xs - S_xs*T'
    # outputs = exodus
  []
  [E0_Va_metal]
    type = ParsedMaterial
    material_property_names = 'T kB'
    constant_names = 'H0_f S0_f'
    constant_expressions = '1.56 3.3'
    function = 'H0_f - S0_f*kB*T' #'-0.2803'#
    f_name = 'E0_Va_metal'
    # outputs = exodus
  []
  [E0_Cr_metal]
    type = ParsedMaterial
    material_property_names = 'T F Gxs'
    function = '(-1572.94 + 157.643*T - 26.908*T*log(T) + 1.89435e-3*T^2 - 1.47721e-6*T^3 + 139250/T '
               '+ Gxs)/F'
    f_name = 'E0_Cr_metal'
    # outputs = exodus
  []

  [E0_Ni_melt]
    type = ParsedMaterial
    f_name = 'E0_Ni_melt'
    constant_names = 'E0_NiF2'
    constant_expressions = '0.473'
    material_property_names = 'E0_Ni_metal kB T n HF_H2 E0_F E_F'
    function = 'E0_Ni_metal + n*(E0_NiF2 + (E_F-E0_F))'
    # outputs = exodus
  []
  [E0_Cr_Cr]
    type = ParsedMaterial
    f_name = 'E0_Cr_Cr'
    material_property_names = 'T F'
    function = '(-8856.94 + 157.48*T - 26.908*T*log(T) + 1.89435e-3*T^2 - 1.47721e-6*T^3 + 139250/T)/F'
  []
  [E0_Cr_melt]
    type = ParsedMaterial
    f_name = 'E0_Cr_melt'
    constant_names = 'E0_CrF2'
    constant_expressions = '-0.39'
    material_property_names = 'E0_Cr_Cr kB T n HF_H2 E0_F E_F'
    function = 'E0_Cr_Cr + n*(E0_CrF2 +  (E_F-E0_F))'
    # outputs = exodus
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
    material_property_names = 'interface_energy_sigma interface_thickness_l'
    function = '3*interface_energy_sigma*interface_thickness_l/4'
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
    args = 'eta1 gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11'# gr12 gr13 gr14'
    material_property_names = 'pi del_int'
    function = '0.5*(1 + tanh( 2*pi*(eta1-0.5)/del_int ) )'
    output_properties = 'h_melt'
    #outputs = exodus
  []
  [switch_metal]
    type = DerivativeParsedMaterial
    f_name = 'h_metal'
    args = 'eta1 gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11'# gr12 gr13 gr14'
    material_property_names = 'pi del_int'
    function = '0.5*(1 - tanh( 2*pi*(eta1 -0.5)/del_int ) )'
    output_properties = 'h_metal'
    #outputs = exodus
  []
  [h_gb]
    type = ParsedMaterial
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 bnds'# gr8 gr9 gr10 gr11 bnds' # gr12 gr13 gr14 bnds'
    f_name = 'h_gb'
    material_property_names = 'pi del_int b l'
    function = '0.5-0.5*tanh((bnds-b)/sqrt(2)/l )'
    # function = '16*( gr0^2*(gr1^2+gr2^2+gr3^2+gr4^2+gr5^2+gr6^2+gr7^2) + '
    #            'gr1^2*(gr2^2+gr3^2+gr4^2+gr5^2+gr6^2+gr7^2) + gr2^3*(gr3^2+gr4^2+gr5^2+gr6^2+gr7^2) '
    #            '+ gr3^2*(gr4^2+gr5^2+gr6^2+gr7^2) '
    #            '+gr4^2*(gr5^2+gr6^2+gr7^2)+gr5^2*(gr6^2+gr7^2)+gr6^2*(gr7^2) )'
    # function = 'tanh(2*pi*bnds*(1-bnds)/del_int )'#'gr^2*eta1^2/(gr^2*eta1^2 + gr^2*eta_melt^2 + eta1^2*eta_melt^2)'
    # function = 'tanh(720*( gr0^2*(gr1^2+gr2^2+gr3^2+gr4^2+gr5^2+gr6^2+gr7^2) + '
    #            'gr1^2*(gr2^2+gr3^2+gr4^2+gr5^2+gr6^2+gr7^2) + '
    #            'gr2^2*(gr3^2+gr4^2+gr5^2+gr6^2+gr7^2) + '
    #            'gr3^2*(gr4^2+gr5^2+gr6^2+gr7^2) +gr4^2*(gr5^2+gr6^2+gr7^2) '
    #            '+gr5^2*(gr6^2+gr7^2) + gr6^2*(gr7^2)    ) )' #'gr^2*eta1^2/(gr^2*eta1^2 + gr^2*eta_melt^2 + eta1^2*eta_melt^2)'
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
    #outputs = exodus
    output_properties = 'c_Ni_metal'
  []
  [c_Ni_melt]
    type = DerivativeParsedMaterial
    f_name = "c_Ni_melt"
    args = 'w_Ni w_Cr'
    material_property_names = 'kB T k_melt E0_Ni_melt E0_Va_melt E0_Cr_melt'
    function = 'exp( (w_Ni - (E0_Ni_melt) - k_melt*kB*T)/kB/T/k_melt )'
    #outputs = exodus
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
    #outputs = exodus
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
    #outputs = exodus
    output_properties = 'c_Ni'
    outputs = nemesis
  []
  [c_Cr]
    type = DerivativeParsedMaterial
    f_name = 'c_Cr'
    material_property_names = 'c_Cr_metal c_Cr_melt h_metal h_melt'
    function = 'c_Cr_metal*h_metal + c_Cr_melt*(1-h_metal)'
    #outputs = exodus
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
    #outputs = exodus
  []
  [del_c_Ni]
    type = DerivativeParsedMaterial
    f_name = 'del_c_Ni'
    args = 'w_Ni w_Cr gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11 eta1 '
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
    args = 'w_Ni w_Cr gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11 eta1 '
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
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Ni w_Cr'# gr8 gr9 gr10 gr11 eta1 w_Ni w_Cr'
    material_property_names = 'kB T k_metal E0_Ni_metal E0_Va_metal  E0_Cr_metal'
    function = 'E0_Va_metal - kB*T*k_metal*log( 1 + exp( (w_Ni - (E0_Ni_metal - '
               'E0_Va_metal))/kB/T/k_metal ) + exp( (w_Cr - (E0_Cr_metal - '
               'E0_Va_metal))/kB/T/k_metal ) )'
    #outputs = exodus
    derivative_order = 2
    outputs = nemesis
    output_properties = 'omega_metal'
  []
  [omega_melt]
    type = DerivativeParsedMaterial
    f_name = 'omega_melt'
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Ni w_Cr'# gr8 gr9 gr10 gr11  eta1 w_Ni w_Cr'
    material_property_names = 'kB T k_melt E0_Ni_melt E0_Va_melt E0_Cr_melt'
    function = 'E0_Va_melt - kB*T*k_melt*(exp( (w_Ni - (E0_Ni_melt - E0_Va_melt) - '
               'k_melt*kB*T)/kB/T/k_melt ) + exp( (w_Cr - (E0_Cr_melt - E0_Va_melt) - '
               'k_melt*kB*T)/kB/T/k_melt ) )'
   outputs = nemesis
   output_properties = 'omega_melt'
    #outputs = exodus
    derivative_order = 2
  []
  [omega_chem]
    type = DerivativeParsedMaterial
    f_name = omega_chem
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Ni w_Cr'# gr8 gr9 gr10 gr11 eta1 w_Ni w_Cr'
    material_property_names = 'h_metal omega_metal omega_melt'
    function = '(omega_metal-omega_melt)'
    outputs = nemesis
    output_properties = 'omega_chem'
    #outputs = exodus
  []
  [susceptibility_Ni]
    type = DerivativeParsedMaterial
    f_name = 'chi_Ni'
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Ni w_Cr'# gr8 gr9 gr10 gr11 eta1  w_Ni w_Cr'
    material_property_names = 'c_Ni_metal(w_Ni,w_Cr) c_Ni_melt(w_Ni,w_Cr) h_metal '
                              'chi_Ni_metal:=D[c_Ni_metal,w_Ni] chi_Ni_melt:=D[c_Ni_melt,w_Ni]'
    function = 'chi_Ni_metal*h_metal + chi_Ni_melt*(1-h_metal)'
    outputs = nemesis
    output_properties = 'chi_Ni'
    #outputs = exodus
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
    #outputs = exodus
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
    #outputs = exodus
  []
  [D_Ni]
    type = ParsedMaterial
    f_name = 'D_Ni'
    material_property_names = 'D_Ni_V D_Ni_GB h_metal h_gb  interface_thickness_l GB_width'
    function = '(D_Ni_V +   1.0*(D_Ni_GB-D_Ni_V )*(GB_width/interface_thickness_l)*h_gb)*h_metal'# + 500*h_melt'
    outputs = nemesis
    output_properties = 'D_Ni'
  []
  [mobility_Ni]
    type = DerivativeParsedMaterial
    f_name = M_Ni
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Ni w_Cr'
    material_property_names = 'chi_Ni h_metal h_melt h_gb D_Ni_V D_Ni_GB interface_thickness_l '
                              'GB_width'
    function = '((D_Ni_V +   1.0*(D_Ni_GB-D_Ni_V )*(GB_width/interface_thickness_l)*h_gb)*h_metal + '
               '500*h_melt )*chi_Ni'
    #outputs = exodus
    outputs = nemesis
    output_properties = 'M_Ni'
    derivative_order = 2
  []

  [susceptibility_Cr]
    type = DerivativeParsedMaterial
    f_name = 'chi_Cr'
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Ni w_Cr'
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
    constant_names = 'D0       D1    E0       E1 '
    constant_expressions = '5.3530e3 19.9  7.8886e7 285000'
    function = '2.2e-7*(c_Va/2.254e-7)'#'1.6e-6*(c_Va/2.254e-7)'
    # function = '(exp(-(84702800778885824000*xc + '
    #            '306299812110263875)/(8927089524736*T))*exp((5885692387260437*xc)/1099511627776)*exp(6'
    #            '99261826406127/35184372088832))*(c_Va/2.254e-7)'
     outputs = nemesis
     output_properties = 'D_Cr_V'
    #outputs = exodus
  []
  [D_Cr_GB]
    type = ParsedMaterial
    f_name = 'D_Cr_GB'
    material_property_names = 'T R xc interface_thickness_l GB_width c_Va'
    constant_names = 'D0       D1    E0       E1'
    constant_expressions = '1.2211e4 17   1.9883e8 2.0e5'
    function = '2.2e-7*7.012e8*(c_Va/2.254e-7)'#'1.6e-6*7.012e6*(c_Va/2.254e-7)'
    # function = '(exp((209787039506469*xc)/17179869184)*exp(2393065153853383/140737488355328)*exp(-(85'
    #            '3955374395308288000*xc + '
    #            '841876566076058375)/(35708358098944*T)))*(c_Va/2.254e-7)/GB_width'
     outputs = nemesis
     output_properties = 'D_Cr_GB'
    #outputs = exodus
  []
  [D_Cr]
    type = ParsedMaterial
    f_name = 'D_Cr'
    material_property_names = 'D_Cr_V D_Cr_GB h_metal h_gb interface_thickness_l GB_width'
    function = '(D_Cr_V +   1.0*(D_Cr_GB-D_Cr_V )*(GB_width/interface_thickness_l)*h_gb)*h_metal'# + 500*h_melt'
    outputs = nemesis
    output_properties = 'D_Cr'
  []
  [mobility_Cr]
    type = DerivativeParsedMaterial
    f_name = M_Cr
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Ni w_Cr'# gr8 gr9 gr10 gr11 eta1  w_Ni w_Cr'
    material_property_names = 'chi_Cr h_metal h_melt h_gb interface_thickness_l GB_width D_Cr_V '
                              'D_Cr_GB'
    function = '((D_Cr_V +   1.0*( D_Cr_GB - D_Cr_V )*(GB_width/interface_thickness_l)*h_gb)*h_metal '
               '+ 500*h_melt )*chi_Cr'
     outputs = nemesis
     output_properties = 'M_Cr'
    #outputs = exodus
    derivative_order = 2
  []
  [M_H]
    type = DerivativeParsedMaterial
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Ni w_Cr phi'# gr8 gr9 gr10 gr11 eta1  w_Ni w_Cr'
    f_name = 'M_H'
    material_property_names = 'h_metal kB T k_melt R'
    constant_names = 'D0_H_melt  E0_H_melt'
    constant_expressions = '2.758e5    36e3 '
    # function = '(1e-7*h_metal + D0_H_melt*exp(-E0_H_melt/R/T)*(1-h_metal) )/kB/T/k_melt'
    function = '(1e-3*h_metal + 1000*(1-h_metal) )/kB/T/k_melt'
    derivative_order = 2
  []
  [M_e]
    type = DerivativeParsedMaterial
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Ni w_Cr phi'# gr8 gr9 gr10 gr11 eta1 '
    f_name = 'M_e'
    material_property_names = 'h_metal h_melt kB T k_metal Va'
    function = '( (48709375000)*h_metal + 1e-12*h_melt )/kB/T/k_metal'
    derivative_order = 2
  []
  #Computing phi
  [flux_Ni]
    type = DerivativeParsedMaterial
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Ni w_Cr phi'# gr8 gr9 gr10 gr11 eta1  w_Ni phi'
    f_name = 'chem_flux_Ni'
    material_property_names = 'z_Ni z_F M_Ni M_H h_metal h_melt'
    function = 'h_melt*M_Ni*2'
    derivative_order = 2
    # function = 'M_Ni'
  []
  [flux_Cr]
    type = DerivativeParsedMaterial
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Cr phi'# gr8 gr9 gr10 gr11 eta1 w_Cr phi'
    f_name = 'chem_flux_Cr'
    material_property_names = 'z_Ni z_F M_Cr M_H h_metal h_melt'
    function = 'h_melt*M_Cr*2'
    #outputs = exodus
    output_properties = 'chem_flux_Cr'
    derivative_order = 2
    # function = 'M_Ni'
  []
  [grad_phi_coeff]
    type = DerivativeParsedMaterial
    f_name = 'elec_flux'
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Cr phi'# gr8 gr9 gr10 gr11 eta1  w_Ni phi'
    material_property_names = 'z_Ni z_F M_Ni M_Cr M_H M_e h_metal h_melt'
    function = 'h_melt*(4*M_Ni + 4*M_Cr + M_H) + (h_metal)*M_e'
    output_properties = 'elec_flux'
    #outputs = exodus
    derivative_order = 2
  []

  #electro_diffusion
  [elec_M_Ni]
    type = DerivativeParsedMaterial
    f_name = 'elec_M_Ni'
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Ni w_Cr phi'# gr8 gr9 gr10 gr11 eta1  phi'
    material_property_names = 'h_metal h_melt M_Ni'
    function = 'h_melt*M_Ni*2'
    derivative_order = 2
  []
  [elec_M_Cr]
    type = DerivativeParsedMaterial
    f_name = 'elec_M_Cr'
    args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 w_Ni w_Cr phi'# gr8 gr9 gr10 gr11 eta1  phi'
    material_property_names = 'h_metal h_melt M_Cr'
    function = 'h_melt*M_Cr*2'
    derivative_order = 2
  []

  #Postprocessor materials
  [metal_thickness]
    type = ParsedMaterial
    f_name = 'metal_thickness'
    # args = 'gr0 gr1'
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

  [elec_conductivity]
    type = ParsedMaterial
    f_name = 'sigma_e'
    # args = 'E_x'
    material_property_names = 'M_e'
    function = '1.43e7*1e6*1e3/1e4' #1.43e7 S/m, *1e6 because electric field is in um, 1e3 to convert A->mA, 1e4 converts /m2 to /cm2
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
  # # [total_void]
  # #   type = ElementIntegralMaterialProperty
  # #   mat_prop = 'is_void'
  # # []
  # [phi_left]
  #   type = PointValue
  #   variable = phi
  #   point = '0 500 0'
  # []
  # [./Jx_wire]
  #   type = FunctionValuePostprocessor
  #   function = 'Jx'
  #   point = '500 950 0'
  # [../]
  # [./Ex_wire]
  #   type = PointValue
  #   point = '500 950 0'
  #   variable = 'E_x'
  # [../]

[]
# [Preconditioning]
#   [SMP]
#     type = SMP
#     full = true
#   []
#   # [No_PC]
#   #   type =
#   # []
#
# []
# [Adaptivity]
#   marker = 'thresh_eta'
#   max_h_level = 3
#   initial_steps = 3
#   # initial_adaptivity = 4
#   [Markers]
#     [thresh_eta]
#       type = ValueThresholdMarker
#       coarsen = 0.95
#       refine = 0.5
#       variable = eta1
#       invert = True
#     []
#   []
# []


[Executioner]
  type = Transient
  solve_type = PJFNK #
  scheme = bdf2
  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm      31                  preonly       lu           2'
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu superlu_dist'
  l_max_its = 30
  l_tol = 1e-3 #1e-10
  # nl_rel_tol = 1e-8
  # nl_abs_tol = 1e-20
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
  dtmin = 1e-7
  end_time = 3.6e6 #1200#8.6e4#20000.0
  automatic_scaling = true
  # compute_scaling_once = false
  # [Adaptivity]
  #   max_h_level = 3#10
  #   coarsen_fraction = 0.05
  #   # initial_adaptivity = 3#10
  #   refine_fraction = 0.9
  #   # weight_names = 'w_Cr w_Ni gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7  eta1'
  #   # weight_values = '1    1   1   1   1   1   1   1   1   1    0.01'
  #   # weight_names = 'w_Cr w_Ni gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7  eta1'
  #   # weight_values = '0.1  0.1   1   1       1     1     1     1     1     1     1'
  # []
  dtmax = 5e3 #500.0
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    growth_factor = 1.25
    cutback_factor = 0.8
  []
  # num_steps = 10
[]

# [Debug]
#   show_var_residual_norms = true
# []
[Outputs]
  # exodus = true
  [nemesis]
    type = Nemesis
    execute_on = 'initial timestep_end final'
    interval = 10
  []
  # nemesis = true
  csv = true
  perf_graph = true
  file_base = 'ni5cr/ni5cr'
[]
