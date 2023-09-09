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
import numpy as np
from string import Template
import re

class ParseEBSD:
    def __init__(self,file_name):
        self.file_name = file_name
        self.extract_header_params()
    def get_header_lines(self):
        with open(self.file_name,'r') as raw_file:
            header_list = []
            while True:
                line = next(raw_file)
                if not (line[0] == '#'):
                    break
                else:
                    header_list+=[line]
        return (header_list)
    def get_points_lines(self):
        with open(self.file_name,'r') as raw_file:
            points_list = []
            while True:
                try:
                    line = next(raw_file)
                    if not (line[0] == '#'):
                        points_list+=[line]
                except:
                    break
        return (points_list)
    def extract_header_params(self):
        header_list = self.get_header_lines()
        self.header_dict = {}
        for line in header_list:
            F = re.match(r"# (.*): (.*)",line)
            if F:
                key,value=F.groups()
                self.header_dict[key]= value
        return self.header_dict
    def is_stackable(self,ebsd_object,axis):
        if axis == "X":
            checklist = ["Y_DIM","Z_DIM","Y_MIN","Y_MAX","Z_MIN","Z_MAX","Y_STEP","Z_STEP","X_STEP"]
            checks = [(self.header_dict[key] == ebsd_object.header_dict[key]) for key in checklist ]
            if all(checks):
                return 1
            else:
                return 0
        else:
            print ("Only 'X' axis stacking has been implemented currently")
            return -1
    def stack_ebsd(self,ebsd_object,axis,file_name):
        if not axis == "X":
            print("Only 'X' axis stacking has been implemented currently")
            return -1
        elif self.is_stackable(ebsd_object,axis):
            #Writing appended ebsd points
            x_bias = float(self.header_dict['X_MAX']) #Adding XMAX through points lines
            grain_num_bias = float(self.header_dict['Features_1'])

            self.header_dict['X_DIM']= str(int(self.header_dict['X_DIM'])+int(ebsd_object.header_dict['X_DIM']))
            self.header_dict['X_MAX'] = str(float(self.header_dict['X_MAX']) +float(ebsd_object.header_dict['X_MAX']))
            self.header_dict['Num_Features'] = str(float(self.header_dict['Num_Features']) +float(ebsd_object.header_dict['Num_Features']))
            with open(file_name,'w+') as new_file:
                #Write header
                for key in self.header_dict.keys():
                    new_file.write("# "+str(key)+ ": "+str(self.header_dict[key]) + "\n#\n" )
                new_file.write("# phi1 PHI phi2 x y z FeatureId PhaseId Symmetry\n")
                #Write self points
                points_lines = self.get_points_lines()
                new_file.writelines(points_lines)


                points_lines = ebsd_object.get_points_lines()

                #Convert to float np array so we can change X values of points
                point_array = np.asarray([line.split() for line in points_lines],np.float32)
                point_array[:,3]+=x_bias #Shift x values
                point_array[:,6]+=grain_num_bias

                #Convert array back into string
                points_lines = "\n".join([" ".join([str(point).rstrip('0').rstrip('.') for point in row]) for row in point_array])
                new_file.writelines(points_lines)

def append_salt(file_name,salt_dim,salt_thickness):

    bias = 0.0
    ebsd_obj = ParseEBSD(file_name)
    H_dict = ebsd_obj.header_dict
    points = ebsd_obj.get_points_lines()

    H_dict[salt_dim+"_MAX"] = str(float(H_dict[salt_dim+"_MAX"]) + salt_thickness)
    N_dim_old = int(H_dict[salt_dim+"_DIM"])
    N_dim_new = math.floor(( float(H_dict[salt_dim+"_MAX"]) - float(H_dict[salt_dim+"_MIN"]) )/float(H_dict[salt_dim+"_STEP"]))
    H_dict[salt_dim+"_DIM"] = str( N_dim_new )

    if (salt_dim == "Z"):
        X_vals = float(H_dict["X_MIN"]) + float(H_dict["X_STEP"])*np.arange(0,int(H_dict["X_DIM"]) )
        Y_vals = float(H_dict["Y_MIN"]) + float(H_dict["Y_STEP"])*np.arange(0,int(H_dict["Y_DIM"]) ) - bias
        Z_vals = float(H_dict["Z_MIN"]) + float(H_dict["Z_STEP"])*np.arange(N_dim_old,int(H_dict["Z_DIM"]) ) - bias
        H_dict["X_MIN"] = str(float(X_vals[0] - bias) )
        H_dict["Y_MIN"] = str(float(Y_vals[0] - bias) )
        H_dict["Z_MIN"] = str(float(H_dict["Z_MIN"]) -bias )
        H_dict["X_MAX"] = str(float(X_vals[-1] + bias) )
        H_dict["Y_MAX"] = str(float(Y_vals[-1] + bias) )
        H_dict["Z_MAX"] = str(float(Z_vals[-1] + bias) )
    if (salt_dim == "X"):
        X_vals = float(H_dict["X_MIN"]) + float(H_dict["X_STEP"])*np.arange(N_dim_old,int(H_dict["X_DIM"]) )
        Y_vals = float(H_dict["Y_MIN"]) + float(H_dict["Y_STEP"])*np.arange(0,int(H_dict["Y_DIM"]) ) - bias
        Z_vals = float(H_dict["Z_MIN"]) + float(H_dict["Z_STEP"])*np.arange(0,int(H_dict["Z_DIM"]) ) - bias
        H_dict["X_MIN"] = str(float(H_dict["X_MIN"]) - bias )
        H_dict["Y_MIN"] = str(float(Y_vals[0] - bias) )
        H_dict["Z_MIN"] = str(float(Z_vals[0] - bias) )
        H_dict["X_MAX"] = str(float(X_vals[-1] + bias) )
        H_dict["Y_MAX"] = str(float(Y_vals[-1] + bias) )
        H_dict["Z_MAX"] = str(float(Z_vals[-1] + bias) )

        #Remove if not running 2D
        H_dict["Z_DIM"] = 0

    # print(min(X_vals),max(X_vals),min(Y_vals),max(Y_vals),min(Z_vals),max(Z_vals))
    G = np.vstack(np.meshgrid(X_vals,Y_vals,Z_vals)).reshape(3,-1).T
    # print(G.shape, 50*50*20)

    S=""

    ###Add new header
    for key in H_dict.keys():
        S = S +("# "+str(key)+ ": "+str(H_dict[key]) + "\n#\n" )

    # print(S)
    ####Append old points
    S = S + "".join(points) + "\n"

    ####Append new points
    P_id = str(0)
    symmetry = H_dict["Symmetry_1"]
    feature_id = str(0)
    phi1 = str(0)
    phi2 = str(0)
    phi3 = str(0)

    # print(G.shape)
    for i in range(G.shape[0]):
        x,y,z = G[i]
        line = " ".join([phi1,phi2,phi3,str(x),str(y),str(z),feature_id,P_id,symmetry] )
        S+=line+"\n"

    with open(file_name,'w') as write_file:
        write_file.write(S)

def single_cubic_phase_equiaxed(xmax,ymax,zmax,h,E,SD,filename):
    nx = int(xmax/h)
    ny = int(ymax/h)
    nz = int(2*E/h)
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
                                         0, 0, 0, nx-1, ny-1, 0, False, False, True,
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
    single_cubic_phase_equiaxed(xmax=150,ymax=30,zmax=2*%G_alloy%,h=0.5,E=%G_alloy%,SD=0.1*%G_alloy%,filename='alloy_grains')
    #Generate alloy grains
    single_cubic_phase_equiaxed(xmax=%coating_thickness%,ymax=30,zmax=2*%G_coating%,h=0.5,E=%G_coating%,SD=0.1*%G_coating%,filename='coating_grains')
    ebsd_1 = ParseEBSD('alloy_grains.txt')
    ebsd_2 = ParseEBSD('coating_grains.txt')

    ebsd_1.stack_ebsd(ebsd_2,'X','EBSD_IC.txt')

    del ebsd_1
    del ebsd_2

    append_salt('EBSD_IC.txt','X',10.0)
