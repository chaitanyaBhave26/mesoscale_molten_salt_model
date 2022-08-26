## Simulates GB diffusion along a GB with 2 grains. By comparing model with Fisher predictions we can understand the effect of interface width and h_GB function on accuracy

[Mesh] #Sets mesh size to 10 microns by 10 microns
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${fparse 150 / 2  }
    ny = ${fparse 30 / 2  }
    xmax = 150
    ymax = 30
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 3
[]

[GlobalParams]
  int_width = 2
  enable_jit = false
  enable_ad_cache = false
[]

[Variables]
  [./c]   #concentration variable
  [../]
[]

[AuxVariables] #Creates 2 grains to interpolate diffusivity using
  [./grain1]
  [../]
  [./grain2]
  [../]
[]

[ICs] #Sets the IC for the second constant phase
  [grain1_IC]
    type = BoundingBoxIC
    variable = grain1
    inside = 1
    outside = 0
    x1 = -10
    y1 = -10
    x2 = 160
    y2 = 15
  []
  [grain2_IC]
    type = BoundingBoxIC
    variable = grain2
    inside = 1
    outside = 0
    x1 = -10
    y1 = 15
    x2 = 160
    y2 = 40
  []
  [c_IC]
    type = FunctionIC
    variable = 'c'
    function = '0.0561'
  []
[]


[Kernels]
  [dt_c]
    type = TimeDerivative
  variable = c
  []
  [./cres]
    type = MatDiffusion
    variable = c
    diffusivity = D
  [../]
[]

[BCs]
  [./right_c] #Fix temperature on the left side
    type = FunctionDirichletBC
    variable = c
    boundary = right
    function = '2.5e-5'
  [../]
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'gamma    gr_energy_sigma  interface_energy_sigma    interface_thickness_l  Va pi '
                 'del_int Na xc GB_width'
    prop_values = '1.5     10                   5.1e-5                        2           1.1131e-11 '
                  ' 3.141592653 ${fparse 0.025 * 2 / 2 } 6.02214076e23 0.0 5e-4'
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
    material_property_names = 'T R xc c_Va'
    function = '(exp(-(84702800778885824000*xc + 306299812110263875)/(8927089524736*T))*exp((5885692387260437*xc)/1099511627776)*exp(699261826406127/35184372088832) )'#*(c_Va/2.254e-7)'
    outputs = exodus
  []
  [D_Cr_GB]
    type = DerivativeParsedMaterial
    f_name = 'D_Cr_GB'
    args = 'grain1 c'
    material_property_names = 'T R xc interface_thickness_l GB_width c_Va'
    function = '(exp((209787039506469*xc)/17179869184)*exp(2393065153853383/140737488355328)*exp(-(853955374395308288000*xc + 841876566076058375)/(35708358098944*T)))'
    outputs = exodus
  []
  [./diffusivity]
    type = DerivativeParsedMaterial
    args = 'grain1 c'
    material_property_names = 'h_melt h_gb D_Cr_V D_Cr_GB GB_width interface_thickness_l'
    function = '(1-h_melt)*(D_Cr_V + h_gb*(  D_Cr_GB/GB_width - D_Cr_V)*(GB_width/interface_thickness_l) ) + h_melt*500'
    outputs = exodus
    f_name = D
  [../]
  [./h_melt]
    type = DerivativeParsedMaterial
    args = 'grain1 grain2'
    f_name = h_melt
    material_property_names = 'pi del_int'
    function = '0' #'eta1:=1-grain1-grain2;0.5*(1+tanh(2*pi*(eta1-0.5)/del_int ))'
    outputs = exodus
  [../]
  [./h_gb]
    type = DerivativeParsedMaterial
    args = 'grain1 grain2'
    f_name = h_gb
    function = '16*grain1^2*grain2^2'
    outputs = exodus
  [../]
[]

[Postprocessors]
  [./right_c]
    type = SideAverageValue
    variable = c
    boundary = right
  [../]
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

[VectorPostprocessors]
  # avoid sampling an element variable on faces
  [./along_gb]
    type = LineValueSampler
    variable = 'c'
    start_point = '0.00 25.0 0'
    end_point = '150 25.0 0'
    num_points = ${fparse 160 / 2 }
    sort_by = id
  [../]
  [./perp_gb]
    type = LineValueSampler
    variable = 'c'
    start_point = '125.00 0.0 0'
    end_point = '125 30.0 0'
    num_points = ${fparse 30 / 2 }
    sort_by = id
  [../]
[]
[Executioner]
  type = Transient

  automatic_scaling = true
  l_max_its = 150
  solve_type = PJFNK#NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'
  l_tol = 1e-04
  nl_rel_tol = 1e-08
  nl_abs_tol = 1e-12
  end_time = 1e5
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
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
