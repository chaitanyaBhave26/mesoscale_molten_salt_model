[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = 'EBSD_IC.txt'
    # filename = '/home/bhavcv/projects/yellowjacket_paper/simulations/ebsd_corr/Ni20Cr_mesh/EBSD_IC.txt' #'test.txt'#
  []
  uniform_refine = 2
  # uniform_refine = 4#2
  parallel_type = DISTRIBUTED
[]

[GlobalParams]
  op_num = 8
  var_name_base = gr
  enable_jit = false
  enable_ad_cache = false
[]

[UserObjects]
  [ebsd_reader]
    type = EBSDReader
    force_preic = true
    execute_on = INITIAL
  []
  [ebsd]
    type = PolycrystalEBSD
    coloring_algorithm = bt
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
    phase = '0'
    execute_on = INITIAL
    # threshold = 0.5
    # coloring_algorithm = jp
  []
  [molten_salt]
    type = PolycrystalEBSD
    coloring_algorithm = jp#bt
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
    phase = '1'

    variable = eta1
    # threshold = 0.5
    # coloring_algorithm = jp
  []
  # [read_exodus]
  #   type = SolutionUserObject
  #   mesh = 'test.e'
  #   timestep = FINAL
  #
  # []
  [grain_tracker]
    type = GrainTracker
    flood_entity_type = ELEMENTAL
    compute_halo_maps = true # For displaying HALO fields
    polycrystal_ic_uo = ebsd
    halo_level = 2 #6
    # threshold = 1e-3
    # threshold = 0.1
    # connecting_threshold = 0.09
    # halo_level = 6
  []
  # [end_sim]
  #   type = Terminator
  #   expression = 'grain_tracker<=469'
  # []
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
    polycrystal_ic_uo = molten_salt
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
  [unique_grains_ic]
    order = CONSTANT
    family = MONOMIAL
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [ghost_elements]
    order = CONSTANT
    family = MONOMIAL
  []
  [halos]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices_ic]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
  [ebsd_grains]
    family = MONOMIAL
    order = CONSTANT
  []
[]

# [Kernels]
#
#
# []

[AuxKernels]
  [BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
  # [ghost_elements]
  #   type = FeatureFloodCountAux
  #   variable = ghost_elements
  #   field_display = GHOSTED_ENTITIES
  #   execute_on = 'initial timestep_end'
  #   flood_counter = grain_tracker
  # []
  [halos]
    type = FeatureFloodCountAux
    variable = halos
    field_display = HALOS
    execute_on = 'initial timestep_end'
    flood_counter = grain_tracker
  []
  [var_indices_ic]
    type = FeatureFloodCountAux
    variable = var_indices_ic
    execute_on = 'initial'
    flood_counter = ebsd
    field_display = VARIABLE_COLORING
  []
  [unique_grains_ic]
    type = FeatureFloodCountAux
    variable = unique_grains_ic
    execute_on = 'initial'
    flood_counter = ebsd
    field_display = UNIQUE_REGION
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
  [grain_aux]
    type = EBSDReaderPointDataAux
    variable = ebsd_grains
    ebsd_reader = ebsd_reader
    data_name = 'feature_id'
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
    # [EulerAngles2RGB]
    #   crystal_structure = cubic
    #   euler_angle_provider = ebsd_reader
    #   grain_tracker = grain_tracker
    # []
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'gamma    gr_energy_sigma  interface_energy_sigma    interface_thickness_l  Va pi '
                 'del_int'
    prop_values = '1.5     6.803                  14.29                        0.5         '
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
    # function = '((del_int*i0*(m/D)^2)/(pi*n*F*R*T*interface_thickness_l*1e-6) )*(m3_to_um3/J_to_eV)/Va' #/Va to normalize L for grand potential in eV/atom rather than eV/um3
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
    all_etas = 'eta1 gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11 '
    phase_etas = 'eta1'
    # outputs = exodus
    # output_properties = 'h_melt'
  []
  [switch_metal]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = h_metal
    all_etas = 'eta1 gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11'
    phase_etas = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11 '
    # outputs = exodus
    # output_properties = 'h_metal'
  []
  # [switch_melt]
  #   type = DerivativeParsedMaterial
  #   f_name = 'h_melt'
  #   args = 'eta1 '
  # []
  # [./switch_void]
  #   type = SwitchingFunctionMultiPhaseMaterial
  #   h_name = h_void
  #   all_etas = 'eta1 gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
  #   phase_etas = 'eta2'
  # [../]
  # [h_gb]
  #   type = ParsedMaterial
  #   args = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 eta1 bnds'
  #   f_name = 'h_gb'
  #   material_property_names = 'pi del_int'
  #   function = 'tanh(2*pi*bnds*(1-bnds)/del_int )' #'gr^2*eta1^2/(gr^2*eta1^2 + gr^2*eta_melt^2 + eta1^2*eta_melt^2)'
  #   # function = 'tanh(720*( gr0^2*(gr1^2+gr2^2+gr3^2+gr4^2+gr5^2+gr6^2+gr7^2) + '
  #   #            'gr1^2*(gr2^2+gr3^2+gr4^2+gr5^2+gr6^2+gr7^2) + '
  #   #            'gr2^2*(gr3^2+gr4^2+gr5^2+gr6^2+gr7^2) + '
  #   #            'gr3^2*(gr4^2+gr5^2+gr6^2+gr7^2) +gr4^2*(gr5^2+gr6^2+gr7^2) '
  #   #            '+gr5^2*(gr6^2+gr7^2) + gr6^2*(gr7^2)    ) )' #'gr^2*eta1^2/(gr^2*eta1^2 + gr^2*eta_melt^2 + eta1^2*eta_melt^2)'
  #   outputs = exodus
  #   output_properties = 'h_gb'
  # []

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
    args = 'w_Ni eta1 gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11' # gr8 gr9'# gr10 gr11 gr12 gr13 gr14'
    material_property_names = 'E0_Ni_metal E0_Ni_melt h_metal'
    function = 'h_metal/E0_Ni_metal + (1-h_metal)/E0_Ni_melt'
  []
  [M_Ni]
    type = DerivativeParsedMaterial
    f_name = 'M_Ni'
    args = 'w_Ni eta1 gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11' # gr8 gr9'# gr10 gr11 gr12 gr13 gr14'
    material_property_names = 'chi_Ni'
    function = '1e-4*chi_Ni'
  []
  [omega_metal]
    type = DerivativeParsedMaterial
    f_name = 'omega_metal'
    args = 'w_Ni eta1 gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11' # gr8 gr9'# gr10 gr11 gr12 gr13 gr14'
    material_property_names = 'E0_Ni_metal c0_Ni_metal'
    function = '-0.5*w_Ni*(w_Ni/E0_Ni_metal + c0_Ni_metal)'
  []
  [omega_melt]
    type = DerivativeParsedMaterial
    f_name = 'omega_melt'
    args = 'w_Ni eta1 gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'# gr8 gr9 gr10 gr11'
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
  # petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold'
  # petsc_options_value = 'hypre    boomeramg      0.7'

  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm      31                  preonly       lu           2'
  l_tol = 1.0e-4
  l_max_its = 20
  nl_max_its = 20
  nl_rel_tol = 1.0e-8

  start_time = 0.0
  end_time = 500
  # num_steps = 30

  [TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.9
    dt = 1.0
    growth_factor = 1.1
    optimal_iterations = 7
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

# [VectorPostprocessors]
#   [mem]
#     type = VectorMemoryUsage
#     execute_on = 'INITIAL TIMESTEP_END NONLINEAR LINEAR'
#     report_peak_value = true
#     mem_units = gigabytes # or bytes, megabytes, gigabytes
#   []
# []
[Outputs]
  # nemesis = true
  [exodus]
    type = Exodus
    execute_on = 'initial timestep_end final'
    interval = 20
  []
  file_base = 'ebsd_reader/ebsd_reader'
  # checkpoint = true
  perf_graph = true
  csv = true
[]
