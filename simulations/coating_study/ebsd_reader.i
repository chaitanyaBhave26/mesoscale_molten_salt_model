[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = 'EBSD_IC.txt'
  []
  uniform_refine = 2
  parallel_type = DISTRIBUTED
[]
# [MultiApps]
#   [corr_simulation]
#     type = FullSolveMultiApp
#     execute_on = FINAL
#     positions = '0 0 0'
#     input_files = 'ni20cr_corr.i'
#     # cli_args = 'GlobalParams/coating_thickness=${GlobalParams/coating_thickness}'
#   []
# []
#
# [Transfers]
#   [to_corr_sim]
#     source_variable = '${GlobalParams/poly_grains} eta1'
#     variable =        '${GlobalParams/poly_grains} eta1'
#     type = MultiAppMeshFunctionTransfer
#     to_multi_app = corr_simulation
#   []
# []
[GlobalParams]
  op_num = 8
  var_name_base = gr
  poly_grains = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
  enable_jit = false
  enable_ad_cache = false
[]

[UserObjects]
  [ebsd_reader]
    type = EBSDReader
  []
  [ebsd]
    type = PolycrystalEBSD
    coloring_algorithm = bt
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
    phase = 1.0
  []
  [ebsd_salt]
    type = PolycrystalEBSD
    coloring_algorithm = bt
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
    phase = 0
    variable = eta1
  []
  [grain_tracker]
    type = GrainTracker
    flood_entity_type = ELEMENTAL
    compute_halo_maps = true # For displaying HALO fields
    polycrystal_ic_uo = ebsd
    halo_level = 6
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = ebsd
    []
  []
  [eta1]
    type = PolycrystalColoringIC
    variable = eta1
    op_index = 0
    polycrystal_ic_uo = ebsd_salt
  []
[]

[Variables]
  [PolycrystalVariables]
  []
  [eta1]
  []
  [w_Ni]
  []
[]

[AuxVariables]
  [bnds]
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
[]



[AuxKernels]
  [BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    execute_on = 'initial timestep_end'
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    execute_on = 'initial timestep_end'
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
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
    prop_names = 'gamma    gr_energy_sigma  interface_energy_sigma    interface_thickness_l  Va pi '
                 'del_int'
    prop_values = '1.5     6.803                  14.29                       0.5         '
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
    prop_values = '0.9        100         0.1         100'
  []
  [L]
    type = ParsedMaterial
    f_name = 'L'
    material_property_names = 'n F kB T Va interface_thickness_l pi del_int'
    constant_names = 'i0   m          D    J_to_eV    m3_to_um3 R'
    constant_expressions = '5.7 58.693e-3  8900 6.242e+18  1e18      8.314' #SI units
    function = '1e-4'
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
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = h_melt
    all_etas = 'eta1 ${GlobalParams/poly_grains}'# gr8 gr9 gr10 gr11 '
    phase_etas = 'eta1'
  []
  [switch_metal]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = h_metal
    all_etas = 'eta1  ${GlobalParams/poly_grains}'# gr8 gr9 gr10 gr11'
    phase_etas = ' ${GlobalParams/poly_grains}'# gr8 gr9 gr10 gr11 '
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
    args = 'w_Ni eta1  ${GlobalParams/poly_grains}'# gr8 gr9 gr10 gr11' # gr8 gr9'# gr10 gr11 gr12 gr13 gr14'
    material_property_names = 'E0_Ni_metal E0_Ni_melt h_metal'
    function = 'h_metal/E0_Ni_metal + (1-h_metal)/E0_Ni_melt'
  []
  [M_Ni]
    type = DerivativeParsedMaterial
    f_name = 'M_Ni'
    args = 'w_Ni eta1  ${GlobalParams/poly_grains}'# gr8 gr9 gr10 gr11' # gr8 gr9'# gr10 gr11 gr12 gr13 gr14'
    material_property_names = 'chi_Ni'
    function = '1e-4*chi_Ni'
  []
  [omega_metal]
    type = DerivativeParsedMaterial
    f_name = 'omega_metal'
    args = 'w_Ni eta1  ${GlobalParams/poly_grains}'# gr8 gr9 gr10 gr11' # gr8 gr9'# gr10 gr11 gr12 gr13 gr14'
    material_property_names = 'E0_Ni_metal c0_Ni_metal'
    function = '-0.5*w_Ni*(w_Ni/E0_Ni_metal + c0_Ni_metal)'
  []
  [omega_melt]
    type = DerivativeParsedMaterial
    f_name = 'omega_melt'
    args = 'w_Ni eta1  ${GlobalParams/poly_grains}'# gr8 gr9 gr10 gr11'
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
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK

  # petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold'
  # petsc_options_value = 'hypre    boomeramg      0.7'
  petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl  -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor'
  petsc_options_value = 'hypre    boomeramg      0.7                                   4                          5  25  HMIS  ext+i  2  0.3'

  l_tol = 1.0e-4
  l_max_its = 20
  nl_max_its = 10
  nl_rel_tol = 1.0e-8

  start_time = 0.0
  end_time = 100
  # num_steps = 1

  [TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.9
    dt = 1.0
    growth_factor = 1.1
    optimal_iterations = 7
  []

  # [Adaptivity]
  #   initial_adaptivity = 2
  #   refine_fraction = 0.7
  #   coarsen_fraction = 0.1
  #   max_h_level = 2
  # []

[]

[Outputs]
  file_base = 'ebsd_reader/ebsd_reader'
  nemesis = true
  csv = true
  # checkpoint = true
  perf_graph = true
[]
