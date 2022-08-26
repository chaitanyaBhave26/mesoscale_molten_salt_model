[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = ${fparse ${xmax} * 3 / ${GlobalParams/int_width} } #250
  xmax = 2.5
  ny = ${fparse ${ymax} * 3 / ${GlobalParams/int_width} }#200
  ymax = 2
  # uniform_refine = 3
[]

[GlobalParams]
  # profile = TANH
  enable_jit = false
  enable_ad_cache = false
  int_width = 0.05
  op_num = 2
  var_name_base = gr
[]

[Variables]
  [w_Ni]
  []

  [gr0] #represents the metal phase
  []
  [gr1]
  []
  [eta1]
  []

[]

[ICs]
  [gr0_IC]
    type = BoundingBoxIC
    variable = 'gr0'
    x1 = -10
    y1 = -10
    x2 = 2.0
    y2 = 1.0
    inside = 1.0
    outside = 0.0
  []
  [gr1_IC]
    type = BoundingBoxIC
    variable = 'gr1'
    x1 = -10
    y1 = 1.0
    x2 = 2.0
    y2 = 10.0
    inside = 1.0
    outside = 0.0
  []
  [eta1_IC]
    type = BoundingBoxIC
    variable = 'eta1'
    x1 = 2.0
    y1 = -10.0
    x2 = 10
    y2 = 10
    inside = 1.0
    outside = 0.0
  []
[]

[AuxVariables]
  [E_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [bnds]
  []
[]

[AuxKernels]
  [BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
[]

[Modules]
  [PhaseField]
    [GrandPotential]
      switching_function_names = 'h_metal h_melt' #eta0 is metal, eta1 is electrolyte
      anisotropic = 'false' #'false false'

      chemical_potentials = 'w_Ni'
      mobilities = 'M_Ni'
      susceptibilities = 'chi_Ni'
      free_energies_w = 'c_Ni_metal c_Ni_melt'

      gamma_gr = g_gr
      mobility_name_gr = L
      kappa_gr = kappa
      free_energies_gr = 'omega_metal omega_melt'

      additional_ops = 'eta1'
      gamma_grxop = g_s
      mobility_name_op = L
      kappa_op = kappa
      free_energies_op = 'omega_metal omega_melt'

    []
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'gamma    gr_energy_sigma  interface_energy_sigma    interface_thickness_l  Va pi '
                 'del_int'
    prop_values = '1.5     6.803                  10                        ${GlobalParams/int_width}         '
                  '1.0951044671410821e-11  3.141592653 0.1'
  []
  [energy_constants]
    type = GenericConstantMaterial
    prop_names = 'kB            T     n  F             k_metal  k_melt   HF_H2  E0_F   E_F' #E0_F is electrode potential in Baes, E_F is pontential in the experimental salt
    prop_values = '8.6173324e-5  973  2  96485.33212   1.0     1.0      1e-9   2.871  3.3607' #12e-4
  []
  [Conserved_param_constants]
    type = GenericConstantMaterial
    prop_names = 'c0_Ni_metal E0_Ni_metal c0_Ni_melt  E0_Ni_melt'
    prop_values = '0.9        1e4         0.1         1e4'
  []
  [L]
    type = ParsedMaterial
    f_name = 'L'
    material_property_names = 'n F kB T Va interface_thickness_l pi del_int'
    constant_names = 'i0   m          D    J_to_eV    m3_to_um3 R'
    constant_expressions = '5.7 58.693e-3  8900 6.242e+18  1e18      8.314' #SI units
    # function = '((del_int*i0*(m/D)^2)/(pi*n*F*R*T*interface_thickness_l*1e-6) )*(m3_to_um3/J_to_eV)/Va' #/Va to normalize L for grand potential in eV/atom rather than eV/um3
    function = '1e-3'
  []

  [./iface]
   # reproduce the parameters from GrandPotentialMultiphase.i
   type = GrandPotentialInterface
   gamma_names = 'g_gr g_s'
   kappa_name = 'kappa'
   sigma       = '6.803  10' # Ratio of 1:1.307 to obtain dihedral angle of 135deg
   width       = ${GlobalParams/int_width}
 [../]
  # #PARAMETERS
  # [kappa] #assume that three interfaces having the same interfacial energy and thickness
  #   type = ParsedMaterial
  #   f_name = kappa_s
  #   material_property_names = 'interface_energy_sigma interface_thickness_l'
  #   function = '3*interface_energy_sigma*interface_thickness_l/4'
  # []
  [kappa_gr] #assume that three interfaces having the same interfacial energy and thickness
    type = ParsedMaterial
    f_name = kappa_gr
    material_property_names = 'gr_energy_sigma interface_thickness_l'
    function = '3*gr_energy_sigma*interface_thickness_l/4'
  []
  # [m]
  #   type = ParsedMaterial
  #   f_name = mu
  #   material_property_names = 'interface_energy_sigma interface_thickness_l'
  #   function = '6*interface_energy_sigma/interface_thickness_l'
  # []
  #SWITCHING FUNCTIONS
  [switch_melt]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = h_melt
    all_etas = 'eta1 gr0 gr1'#
    phase_etas = 'eta1'
    # outputs = exodus
    # output_properties = 'h_melt'
  []
  [switch_metal]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = h_metal
    all_etas = 'eta1 gr0 gr1'# gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11'
    phase_etas = 'gr0 gr1'# gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11 '
  []

  [c_Ni_metal]
    type = DerivativeParsedMaterial
    f_name = 'c_Ni_metal'
    args = 'w_Ni'
    material_property_names = 'kB T k_metal E0_Ni_metal c0_Ni_metal'
    function = 'w_Ni/E0_Ni_metal + c0_Ni_metal'
  []
  [c_Ni_melt]
    type = DerivativeParsedMaterial
    f_name = 'c_Ni_melt'
    args = 'w_Ni'
    material_property_names = 'kB T k_melt E0_Ni_melt c0_Ni_melt'
    function = 'w_Ni/E0_Ni_melt + c0_Ni_melt'
  []
  [chi_Ni]
    type = DerivativeParsedMaterial
    f_name = 'chi_Ni'
    args = 'w_Ni eta1 gr0 gr1'# gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11' # gr8 gr9'# gr10 gr11 gr12 gr13 gr14'
    material_property_names = 'E0_Ni_metal E0_Ni_melt h_metal'
    function = 'h_metal/E0_Ni_metal + (1-h_metal)/E0_Ni_melt'
  []
  [M_Ni]
    type = DerivativeParsedMaterial
    f_name = 'M_Ni'
    args = 'w_Ni eta1 gr0 gr1'# gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11' # gr8 gr9'# gr10 gr11 gr12 gr13 gr14'
    material_property_names = 'chi_Ni'
    function = '1e-4*chi_Ni'
  []
  [omega_metal]
    type = DerivativeParsedMaterial
    f_name = 'omega_metal'
    args = 'w_Ni eta1 gr0 gr1'# gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11' # gr8 gr9'# gr10 gr11 gr12 gr13 gr14'
    material_property_names = 'E0_Ni_metal c0_Ni_metal'
    function = '-0.5*w_Ni*(w_Ni/E0_Ni_metal + c0_Ni_metal)'
  []
  [omega_melt]
    type = DerivativeParsedMaterial
    f_name = 'omega_melt'
    args = 'w_Ni eta1 gr0 gr1'# gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11'
        material_property_names = 'E0_Ni_melt c0_Ni_melt'
    function = '-0.5*w_Ni*(w_Ni/E0_Ni_melt + c0_Ni_melt)'
  []
  [gb_perimeter]
    type = ParsedMaterial
    f_name ='gb_perimeter'
    args = 'bnds'
    material_property_names = 'interface_thickness_l'
    function = '(bnds>=0.08 & bnds<=0.5)/interface_thickness_l'
  []

[]

[Postprocessors]
  [dt]
    type = TimestepSize
  []
  [n_elements]
    type = NumElems
    execute_on = 'initial timestep_end'
  []
  [n_nodes]
    type = NumNodes
    execute_on = 'initial timestep_end'
  []
  [DOFs]
    type = NumDOFs
  []
  [GB_length]
    type = ElementIntegralMaterialProperty
    mat_prop = 'gb_perimeter'
  []
[]
# [Preconditioning]
#   [SMP]
#     type = SMP
#     full = true
#   []
# []
[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  automatic_scaling = true
  petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre    boomeramg      0.7'

  # petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  # petsc_options_value = 'asm      31                  preonly       lu           2'
  l_tol = 1.0e-4
  l_max_its = 20
  nl_max_its = 20
  nl_abs_tol=1e-8
  nl_rel_tol = 1.0e-12

  start_time = 0.0
  end_time = 1e6
  steady_state_detection = true
  steady_state_tolerance = 1e-12
  # num_steps = 30

  [TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.8
    dt = 0.1
    growth_factor = 1.25
    optimal_iterations = 5
  []

  # [Adaptivity]
  #   initial_adaptivity = 2
  #   refine_fraction = 0.9
  #   coarsen_fraction = 0.05
  #   max_h_level = 2
  # []
  # num_steps = 1
[]
# [Debug]
#   show_var_residual_norms = true
# []


[Outputs]
  # nemesis = true
  [exodus]
    type = Exodus
    execute_on = 'initial timestep_end final'
    interval = 20
  []
  file_base = 'gamma_ic/gamma_ic'
  perf_graph = true
  csv = true
[]
