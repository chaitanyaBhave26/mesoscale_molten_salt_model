## Simulates GB diffusion along a GB with 2 grains. By comparing model with Fisher predictions we can understand the effect of interface width and h_GB function on accuracy

[Mesh] #Sets mesh size to 10 microns by 10 microns
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = '${fparse ${xmax} / %int_width%  }'
    ny = '${fparse ${ymax} / %int_width%  }'
    xmax = 150
    ymax = 12 #30
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 2
[]

[GlobalParams]
  # int_width = %int_width%
  enable_jit = false
  enable_ad_cache = false
[]

[Variables]
  [c] #concentration variable
  []
[]

[AuxVariables] #Creates 2 grains to interpolate diffusivity using
  [grain1]
  []
  [grain2]
  []
[]

[ICs] #Sets the IC for the second constant phase
  [grain1_IC]
    type = FunctionIC
    variable = grain1
    function = 'y0:=${Mesh/gen/ymax} / 2;0.5-0.5*tanh(2*(y-y0)/%int_width%)'
  []
  [grain2_IC]
    type = FunctionIC
    variable = grain2
    function = 'y0:=${Mesh/gen/ymax} / 2;0.5+0.5*tanh(2*(y-y0)/%int_width%)'
  []
  [c_IC]
    type = FunctionIC
    variable = 'c'
    function = '0'
  []
[]

[Kernels]
  [dt_c]
    type = TimeDerivative
    variable = c
  []
  [cres]
    type = MatDiffusion
    variable = c
    diffusivity = D
  []
[]

[BCs]
  [right_c] #Fix temperature on the left side
    type = FunctionDirichletBC
    variable = c
    boundary = right
    function = '1'
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'gamma    gr_energy_sigma  interface_energy_sigma    interface_thickness_l  Va pi '
                 'del_int Na xc GB_width b l'
    prop_values = '1.5     10                   5.1e-5                        %int_width%           '
                  '1.1131e-11  3.141592653 ${fparse 0.025 * 2 / %int_width% } 6.02214076e23 0.0 5e-4 '
                  '0.8 0.01'
  []
  [energy_constants]
    type = GenericConstantMaterial
    prop_names = 'kB            R      T     n  F             k_metal  k_melt   HF_H2  E0_F   E_F' #E0_F is electrode potential in Baes, E_F is pontential in the experimental salt
    prop_values = '8.6173324e-5 8.314  973   2  96485.33212   1.0     1.0      1e-9   2.871  3.3607'
  []

  [D_Cr_V]
    type = DerivativeParsedMaterial
    f_name = 'D_Cr_V'
    args = 'grain1 c'
    material_property_names = 'T R xc'
    function = '2.072e-7'
    outputs = exodus
  []
  [D_Cr_GB]
    type = DerivativeParsedMaterial
    f_name = 'D_Cr_GB'
    args = 'grain1 c'
    material_property_names = 'T R xc interface_thickness_l GB_width'
    function = '(7.26e-4)'
    outputs = exodus
  []
  [diffusivity]
    type = DerivativeParsedMaterial
    args = 'grain1 c'
    material_property_names = 'h_melt h_gb D_Cr_V D_Cr_GB GB_width interface_thickness_l'
    function = '(1-h_melt)*(D_Cr_V + h_gb*(  D_Cr_GB/GB_width - '
               'D_Cr_V)*(GB_width/interface_thickness_l) )'
    outputs = exodus
    f_name = D
  []
  [h_melt]
    type = DerivativeParsedMaterial
    args = 'grain1 grain2'
    f_name = h_melt
    material_property_names = 'pi del_int'
    function = 'eta1:=1-grain1-grain2;0.5*(1+tanh(2*pi*(eta1-0.5)/del_int ))'
    outputs = exodus
  []
  [h_gb]
    type = DerivativeParsedMaterial
    args = 'grain1 grain2'
    f_name = h_gb
    material_property_names = 'b l pi'
    function = 'bnds:=grain1^2+grain2^2;0.5-0.5*tanh(pi*(bnds-b)/l )'
    outputs = exodus
  []
[]

[Postprocessors]
  [right_c]
    type = SideAverageValue
    variable = c
    boundary = right
  []
  [total_c]
    type = ElementIntegralVariablePostprocessor
    variable = c
  []
  [avg_D]
    type = ElementAverageMaterialProperty
    mat_prop = 'D'
  []
  [elapsed]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  []
[]

[Executioner]
  type = Transient

  automatic_scaling = true
  l_max_its = 150
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart '
                        '-pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'
  l_tol = 1e-04
  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-12
  end_time = 1e6 #3.6e6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    iteration_window = 2
    optimal_iterations = 9
    growth_factor = 1.2
    cutback_factor = 0.8
  []
  dtmax = 1e4
[]

[Outputs]
  file_base = '2_grain/2_grain'
  [exodus]
    type = Exodus
    execute_on = 'INITIAL TIMESTEP_END FINAL'
    interval = 20
  []
  csv = true
[]
