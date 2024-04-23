# Ivan Arias
# 2021/02/15

import create_1b_from1a as dq_1b
import os
import glob
import pyart
import tempfile
import shutil
import multiprocessing as mp
from datetime import datetime
import gc

# Function parameter is the path where original file are located

input_path = '/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/networked_radar/Enhanced_Coeff/data/original'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/networked_radar/Enhanced_Coeff/data/case_0230/CSAPR/original'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/networked_radar/Enhanced_Coeff/data/case_0230/original'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/networked_radar/Enhanced_Coeff/data/CSAPR/original'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/networked_radar/Enhanced_Coeff/data/original'
'/net/k2/storage/people/idariash/home/tmp/drops_pat/original'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/networked_radar/Enhanced_Coeff/data/original'
'/net/k2/storage/projects/RELAMPAGO/RMA1/Original_data'
'/net/k2/storage/projects/RELAMPAGO/RMA1/20181214'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/EOL/data/rma/original'
'/net/k2/storage/projects/CSAPR/DOE_b1'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/Dec_14/CSAPR/original'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1a/2019'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/test/data_1a'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1a/2019/01'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/test_RainRate/data_1a'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1a/'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/test/data_1a'

output_path = '/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/networked_radar/Enhanced_Coeff/data/drops'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/networked_radar/Enhanced_Coeff/data/case_0230/CSAPR/drops'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/networked_radar/Enhanced_Coeff/data/case_0230/drops'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/networked_radar/Enhanced_Coeff/data/CSAPR/drops'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/networked_radar/Enhanced_Coeff/data/drops'
'/net/k2/storage/people/idariash/home/tmp/drops_pat/original'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/networked_radar/Enhanced_Coeff/data/drops'
'/net/k2/storage/projects/RELAMPAGO/RMA1/quality_controlled/level_1b'
'/net/k2/storage/projects/RELAMPAGO/RMA1/drops/20181214'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/EOL/data/rma/corrected'
'/net/k2/storage/projects/CSAPR/level_1b.2'
'/net/k2/storage/people/idariash/home/Field_Campaigns/Relampago/Analisis/Dec_14/CSAPR/level_1b/test'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/test/data_1b_whitney'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/test/att_correction_test_1b'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/level_1b.2'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/test_RainRate/data_1b'
'/net/k2/storage/projects/RELAMPAGO/quality_controlled_data/test/ForJacob'

log_name = os.path.join(output_path, 'log_1b_data.log')

log = open(log_name,'a')
log.write('1b data generation log ' + '\n')
log.write(str(datetime.now()) + '\n')
log.close()


def create_1b_from_1a(file):
    
    # file_size = os.path.getsize(file)/1E6
    # if file_size < 40:
    #     return 
    
    tmp_path = tempfile.mkdtemp()
    tmp_path_radx = os.path.join(tmp_path, 'radx')
    tmp_path_drops = os.path.join(tmp_path, 'drops')
    os.mkdir(tmp_path_radx)
    os.mkdir(tmp_path_drops)
    
    try:
    
        dq_1b.RadxConvert_1a(file, tmp_path_radx)
        
        dq_1b.run_drops(tmp_path_radx, tmp_path_drops)
        
        radar = dq_1b.join_files(tmp_path_drops)
        
        radar_1a = pyart.io.read(file)
        
        radar = dq_1b.add_fields_from_1a(radar, radar_1a) # for CHIVO
        
        #radar = dq_1b.add_fields_from_b1(radar, radar_1a) # for CSAPR
        
        #radar = dq_1b.add_fields_from_rma_cfradial(radar, radar_1a) # for RMA
        
        dq_1b.exportFile(radar, file, output_path)
        
        del radar, radar_1a
        
    except:
        
        log = open(log_name,'a')
        
        log.write('Error in file ' + file + '\n')
       
        log.close()
        
    gc.collect()
    
    shutil.rmtree(tmp_path)
    




files = [f for f in glob.glob(input_path +"/**/*.nc", recursive=True)]
files.sort()

for file in files:
 	create_1b_from_1a(file)

# pool = mp.Pool(16)
# pool.map_async(create_1b_from_1a, files)
# pool.close()
# pool.join()
