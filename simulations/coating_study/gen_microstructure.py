'''
Test of Single Cubic Phase Equiaxed Pipeline example in Python
'''

from dream3d import simpl
from dream3d import simplpy
import dream3d.simpl_helpers as sh
# from dream3d import simpl_test_dirs as sd
import dream3d.orientationanalysispy as orientationanalysis
import dream3d.statstoolboxpy as statstoolboxpy
import dream3d.syntheticbuildingpy as syntheticbuilding
from dream3d.samplingpy import crop_image_geometry
import math

def single_cubic_phase_equiaxed(xmax,ymax,zmax,h,E,SD,filename):
    nx = int(xmax/h)-1
    ny = int(ymax/h)-1
    nz = int(2*E/h)-1
    mu = math.log(E) - 0.5*math.log( (SD/E)**2 + 1 )
    sd = math.sqrt(math.log( (SD/E)**2 + 1 ))
    # print(nx,ny,nz)

    # Create Data Container Array
    dca = simpl.DataContainerArray()

    # Stats Generator Filter
    # # Create a StatsDataArray
    stats_data_array = simpl.StatsDataArray(1, 'StatsDataArray', True)
    stats_data_array.fillArrayWithNewStatsData(1, simpl.PhaseType.Primary)
    stats_data_array.getStatsData(0).BinStepSize = 10.0

    # Might not need to set the Min and Max Feature Diameter (should be calculated)
    stats_data_array.getStatsData(0).MinFeatureDiameter = max(h*5,E-3*SD)
    stats_data_array.getStatsData(0).MaxFeatureDiameter = min(ny,E+3*SD)

    stats_data_array.getStatsData(0).PhaseFraction = 1
    stats_data_array.getStatsData(0).BoundaryArea = 0
    stats_data_array.getStatsData(0).FeatureSize_DistType = 1

    # Crystal Structures:
    crystal_structures = simpl.UInt32Array(1, 'CrystalStructures', 1)

    # Phase Types
    phase_types = simpl.UInt32Array(1, 'PhaseTypes', simpl.PhaseType.Primary)

    # Phase Name
    phase_names = simpl.StringDataArray(1, 'PhaseName', True)
    phase_names.setValue(0, 'Primary')

    # TESTING --> FloatArray array

    cDims = simpl.VectorSizeT([1])
    floatArray1 = simpl.FloatArray(1, cDims, 'Average', mu)

    floatArray2 = simpl.FloatArray(1, cDims, 'Standard Deviation', sd)

    stats_data_array.getStatsData(0).FeatureSizeDistribution = [floatArray1, floatArray2]

    # Using the GenerateStatsData and CreateDynamicTableData functions
    euler_dynamic_table_data = sh.CreateDynamicTableData([[0, 0, 0, 0, 0]],
                                                      ['Euler 1', 'Euler 2', 'Euler 3', 'Weight', 'Sigma'], ['0'])
    axis_dynamic_table_data = sh.CreateDynamicTableData([[0, 0, 0, 0, 0]],
                                                     ['Angle(w)', 'Axis (h)', 'Axis (k)', 'Axis (l)', 'Weight (MRD)'],
                                                     ['0'])
    err = syntheticbuilding.generate_primary_stats_data(dca, 'Primary', 0, 1, 0,
                               1, mu, sd, 1, 2, 3*sd/10,
                               True, 'StatsGeneratorDataContainer',
                               'CellEnsembleData', False,
                               simpl.DataArrayPath('', '', ''),
                               euler_dynamic_table_data,
                               axis_dynamic_table_data,
                               euler_dynamic_table_data)

    assert err == 0, f'StatsGeneratorFilter ErrorCondition: {err}'

    # Initialize Synthetic Volume
    err = syntheticbuilding.initialize_synthetic_volume(dca, 'SyntheticVolumeDataContainer', 'CellData',
                                                        'CellEnsembleMatrixName',
                                                        6,
                                                        simpl.IntVec3([nx, ny, nz]),
                                                        simpl.FloatVec3([h, h, h]),
                                                        simpl.FloatVec3([0, 0, 0]),
                                                        simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                            'CellEnsembleData', 'Statistics'),
                                                        simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                            'CellEnsembleData', 'PhaseTypes'),
                                                        simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                            'CellEnsembleData', 'PhaseName'),
                                                        False, 0, 'NOT NEEDED')
    assert err == 0, f'InitializeSyntheticVolume ErrorCondition: {err}'

    # Establish Shape Types
    err = syntheticbuilding.establish_shape_types(dca,
                                                  simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                      'CellEnsembleData', 'PhaseTypes'),
                                                  'ShapeTypes', [simpl.ShapeType.Ellipsoid])
    assert err == 0, f'EstablishShapeTypes ErrorCondition: {err}'

    # Pack Primary Phases
    err = syntheticbuilding.pack_primary_phases(dca,
                                                simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData',
                                                                    ''),
                                                'Grain Data', 'CellEnsembleData', 'FeatureIds', 'Phases',
                                                'Phases', 'NumFeatures',
                                                simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                    'CellEnsembleData', 'Statistics'),
                                                simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                    'CellEnsembleData', 'PhaseTypes'),
                                                simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                    'CellEnsembleData', 'PhaseName'),
                                                simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                    'CellEnsembleData', 'ShapeTypes'),
                                                simpl.DataArrayPath('', '', ''), False,
                                                sh.FeatureGeneration.GenerateFeatures, '', '', False,
                                                False, sh.SaveShapeDescArrays.DoNotSave,
                                                simpl.DataArrayPath('', '', ''),
                                                simpl.DataArrayPath('', '', ''))
    assert err == 0, f'PackPrimaryPhases ErrorCondition: {err}'

    # Find Feature Neighbors
    err = statstoolboxpy.find_neighbors(dca, simpl.DataArrayPath('SyntheticVolumeDataContainer', 'Grain Data', ''),
                                    'SharedSurfaceAreaList', 'NeighborList',
                                    simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData',
                                                        'FeatureIds'),
                                    '', 'NumNeighbors', 'SurfaceFeatures', False, True)
    assert err == 0, f'FindNeighbors ErrorCondition: {err}'

    # Match Crystallography
    err = syntheticbuilding.match_crystallography(dca, simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                           'CellEnsembleData', 'Statistics'),
                                                  simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                      'CellEnsembleData', 'CrystalStructures'),
                                                  simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                      'CellEnsembleData', 'PhaseTypes'),
                                                  simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                      'CellData', 'FeatureIds'),
                                                  simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                      'Grain Data', 'Phases'),
                                                  simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                      'Grain Data', 'SurfaceFeatures'),
                                                  simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                      'Grain Data', 'NeighborList'),
                                                  simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                      'Grain Data', 'SharedSurfaceAreaList'),
                                                  simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                      'CellEnsembleData', 'NumFeatures'),
                                                  'EulerAngles', 'Volumes', 'EulerAngles', 'AvgQuats', 100000)
    assert err == 0, f'MatchCrystallography ErrorCondition: {err}'

    # Generate IPF Colors
    err = orientationanalysis.generate_ipf_colors(dca, simpl.FloatVec3([0, 0, 1]),
                                                  simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                      'CellData', 'Phases'),
                                                  simpl.DataArrayPath('SyntheticVolumeDataContainer',
                                                                      'CellData', 'EulerAngles'),
                                                  simpl.DataArrayPath('StatsGeneratorDataContainer',
                                                                      'CellEnsembleData', 'CrystalStructures'),
                                                  False,
                                                  simpl.DataArrayPath('', '', ''), 'IPFColor')
    assert err == 0, f'GenerateIPFColors ErrorCondition: {err}'

    err = crop_image_geometry(dca,'',
                                         simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData', ''),
                                         simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellEnsembleData', ''),
                                         0, 0, 0, xmax, ymax, 0, False, False, True,
                                         simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData', 'FeatureIds'))
    assert err == 0, f'CropImageGeometry ErrorCondition: {err}'

    # Write to DREAM3D file
    # err = sh.WriteDREAM3DFile(filename+'.dream3d',dca)
    # assert err == 0, f'WriteDREAM3DFile ErrorCondition: {err}'



    #Write INL file
    err = orientationanalysis.inl_writer(dca,filename+'.txt',
                                        simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellEnsembleData', 'PhaseName'),
                                        simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData', 'FeatureIds'),
                                        simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData', 'Phases'),
                                        simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellEnsembleData',
                                                            'CrystalStructures'),
                                        simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellEnsembleData', 'NumFeatures'),
                                        simpl.DataArrayPath('SyntheticVolumeDataContainer', 'CellData', 'EulerAngles')
                                        )#,simpl.DataArrayPath('SyntheticVolumeDataContainer','Phases','MaterialName'))

    # err = orientationanalysis.inl_writer(dca, 'EBSD_IC.inl')
    assert err == 0, f'INL Writer ErrorCondition {err}'

if __name__ == '__main__':
    #Generate alloy grains
    single_cubic_phase_equiaxed(xmax=150,ymax=50,zmax=2*%G_alloy%,h=0.5,E=%G_alloy%,SD=0.1*%G_alloy%,filename='alloy_grains')
    #Generate alloy grains
    single_cubic_phase_equiaxed(xmax=%coating_thickness%,ymax=50,zmax=2*%G_coating%,h=0.5,E=%G_coating%,SD=0.1*%G_coating%,filename='coating_grains')
