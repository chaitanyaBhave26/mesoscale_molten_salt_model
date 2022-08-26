[Mesh]
  [gen]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    xmin = 0
    xmax = 160
    nx = 160

    ymin = 0
    ymax = 30
    ny = 30
    elem_type = QUAD4
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 2
[]


[MultiApps]
  [alloy_grains]
    type = FullSolveMultiApp
    # not setting app_type to use the same app type of master, i.e. MooseTestApp
    execute_on = 'INITIAL'
    positions = '0 0 0'
    input_files = 'voronoi_grains.i'
  []
  # [coating_grains]
  #   type = FullSolveMultiApp
  #   # not setting app_type to use the same app type of master, i.e. MooseTestApp
  #   execute_on = 'INITIAL'
  #   positions = '150 0 0'
  #   input_files = 'coating_grains.i'
  # []
[]

[Transfers]
  [from_alloy_grains]
    source_variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13'
    variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13'
    type = MultiAppMeshFunctionTransfer
    from_multi_app = alloy_grains
  []
  # [from_coating_grains]
  #   source_variable = 'gr0  gr1   gr2   gr3   gr4   gr5   gr6   gr7   gr8   gr9'
  #   variable =        'gr14  gr15  gr16  gr17  gr18  gr19  gr20  gr21 gr22 gr23'
  #   type = MultiAppMeshFunctionTransfer
  #   from_multi_app = coating_grains
  # []
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 14 # Number of order parameters used
  var_name_base = gr # Base name of grains
  enable_jit = false
  enable_ad_cache = false
  poly_grains = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13'# gr14 gr15 gr16 gr17 gr18 gr19 gr20 gr21 gr22 gr23'
[]
[UserObjects]
  [./grain_tracker]
    type = GrainTracker
    threshold = 0.2
    connecting_threshold = 0.08
    compute_halo_maps = true # Only necessary for displaying HALOS
  [../]
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  [./PolycrystalVariables]
  [../]

  [eta1]
  []
  [w_Ni]
  []
[]
[AuxVariables]
  # Dependent variables
  [./bnds]
    # Variable used to visualize the grain boundaries in the simulation
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # [./ghost_regions]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./halos]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
[]
[AuxKernels]
  # AuxKernel block, defining the equations used to calculate the auxvars
  [./bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = 'initial timestep_end'
  [../]
  # [./ghosted_entities]
  #   type = FeatureFloodCountAux
  #   variable = ghost_regions
  #   flood_counter = grain_tracker
  #   field_display = GHOSTED_ENTITIES
  #   execute_on = 'initial timestep_end'
  # [../]
  # [./halos]
  #   type = FeatureFloodCountAux
  #   variable = halos
  #   flood_counter = grain_tracker
  #   field_display = HALOS
  #   execute_on = 'initial timestep_end'
  # [../]
[]
[ICs]
  [eta1_IC]
    type = BoundingBoxIC
    variable = 'eta1'
    x1 = 150
    y1 = -10
    x2 = 300
    y2 = 100
    inside = 1
    outside = 0
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
  [GB_length]
    type = ElementIntegralMaterialProperty
    mat_prop = 'gb_perimeter'
  []
[]


[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre boomeramg 101 ds'

  l_max_its = 30 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 40 # Max number of nonlinear iterations
  nl_rel_tol = 1e-10 # Absolute tolerance for nonlienar solves
  nl_abs_tol = 1e-12
  start_time = 0.0
  end_time = 500

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 25 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 # Time step will adapt to maintain this number of nonlinear iterations
  [../]

  # num_steps = 2
[]

[Outputs]
  file_base = 'vor_ic/vor_ic'
  exodus = true
  csv = true
  perf_graph = true
[]
