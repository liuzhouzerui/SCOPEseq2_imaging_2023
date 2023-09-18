# Example workflow 
## project folder arrangement
project

> bead
>> count00000_channelcy3,cy5,bf_seq0000_cy3.tif_Results.xls \
>> count00000_channelcy3,cy5,bf_seq0000_cy5.tif_Results.xls \
>> ... \
>> count00007_channelcy3,cy5,bf_seq0015_cy5.tif_Results.xls \
>> 192_8mer_seq_reference.csv

> well
>> ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_GFP.tif_Results \
>> ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_TRITC.tif_Results \
>> ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_CY5 QUAD.tif_Results \
>> UL_BR_0.xls \
>> *Results.csv* 

> cell
>>  GFP 
>>> intensity_data

>> TRITC
>>> intensity_data

>> CY5
>>> intensity_data

>> single_cell

> seq

> analysis
>> *well_intensity.pdf* \
>> *probe_intensity.pdf* \
>> *register.pdf* \
>> *sc.image.matrix.txt* \
>> *sc.edrops.matrix.txt* \
>> *sc.edrops.matrix.hist.txt*

> *bead_intensity.obj* \
> *cell_intensity.obj* \
> *cell_bead_link.obj* \
> *image_seq_link.obj* 


## analysis workflow

### step 1, measure well intensity, in ImageJ
- INPUT: a folder contains the following tif files from cell imaging
   > ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_BF.tif \
   > ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_GFP.tif \
   > ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_TRITC.tif \
   > ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_CY5 QUAD.tif

- RUN: a ImageJ macro script
   > OBCV2_wells_cells.txt

- OUTPUT: a 'intensity_data' folder contains the following xls files
   > ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_GFP.tif_Results.xls \
   > ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_TRITC.tif_Results.xls \
   > ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_CY5 QUAD.tif_Results.xls

### step 2, get wells with cells, in Python
- INPUT:
1. project_folder - the location of the project.
2. set_device 
   > n_lane: int - index of device lane \
   > total_lanes: int - number of lanes in the device. 1 for small, medium small device. 5 for large device. \
   > d_th: int - radius of well (pixel)
3. set_well_scan
   > well_channel: list - fluorescent channels used in cell scan

- RUN:
```python
from scopeseq.Bead_Intensity import BeadIntensity # for bead and bead-cell-link
from scopeseq.Cell_Intensity import CellIntensity # for cell and bead-cell-link
from scopeseq.Cell_Bead_Link import CellBeadLink # for bead-cell-link and image_seq_link
from scopeseq.Image_Seq_Link import ImageSeqLink # for image_seq_link
from scopeseq.Seq import CountMatrix # for seq
#import anndata # for seq and schpf
#from scopeseq.ScHPF import ScHPF # for schpf

from scopeseq.utils import pickle_dump, pickle_load # for all
from scopeseq.plot import plt_map
from scopeseq.stats import fit_gaussian, fit_double_gaussian
from scopeseq.clustering import cptt

project_folder  = 'your_project_path'
```

```python
cell_intensity = CellIntensity(project_folder)
cell_intensity.set_device(n_lane=0, total_lanes=1, d_th=40)
cell_intensity.set_well_scan(well_channel=['GFP', 'TRITC', 'CY5 QUAD'])
cell_intensity.generate_well_intensity()
```
- OUTPUT:
1. project/analysis
   > well_intensity.pdf

- INPUT:
1. cell_filter
   > filter_th: dict - {fluorescent channel: threshold}. threshold get from well_intensity.pdf
   
- RUN: 
```python
cell_intensity.cell_filter(filter_th={'TRITC': 2**9.5, 'CY5 QUAD': 2**9.3})
pickle_dump(cell_intensity, project_folder + 'cell_intensity.obj')
```

- OUTPUT:
1. project/well
   > well_lite.csv

### step 3, single cell intensity measurement, in imageJ
- INPUT: 
1. the output from step 2, well_lite.csv
2. a folder contains the following tif files from cell imaging
   > ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_GFP.tif \
   > ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_TRITC.tif \
   > ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_CY5 QUAD.tif

- RUN:
   > OBCV2_well_cell_crop \
   > OBCV2_well_cell_merge \
   > OBCV2_well_cell_measure

- OUTPUT:
1. intensities data folders
   > GFP/intensity_data \
   > TRITC/intensity_data \
   > CY5/intensity_data
2. single cell identification folder - single_cell

- RUN:
   > manually delete doublet images in folder, single_cell.


### step 4, generate single cell fluorescent intensities, in Python
- INPUT: 
1. the output from step 3, intensities data folders
2. the output from step 3, the single_cell folder
3. cell_croped_image_channels - cell image channels

- RUN:
```python
cell_intensity = pickle_load(project_folder + 'cell_intensity.obj')
cell_intensity.generate_cell_intensity(cell_channel=['GFP', 'TRITC', 'CY5'])
pickle_dump(cell_intensity, project_folder + 'cell_intensity.obj')
```

- OUTPUT:
1. project/cell
   > sc_image.matrix.txt

### step 5, bead intensity measurement, in imageJ
- INPUT: a folder contains the following tif files from bead imaging
   > count00000_channelcy3,cy5,bf_seq0000_bf.tif \
   > count00000_channelcy3,cy5,bf_seq0000_cy3.tif \
   > count00000_channelcy3,cy5,bf_seq0000_cy5.tif \
   > ... \
   > count00007_channelcy3,cy5,bf_seq0015_cy5.tif 

- RUN: a ImageJ macro script
   > OBCV2_bead_identification_intensity_measure.txt

- OUTPUT: a 'intensity_data' folder contains the following xls files 
   > count00000_channelcy3,cy5,bf_seq0000_cy3.tif_Results.xls \
   > count00000_channelcy3,cy5,bf_seq0000_cy5.tif_Results.xls \
   > ... \
   > count00007_channelcy3,cy5,bf_seq0015_cy5.tif_Results.xls 

- NOTE: the filename format 
   > 00000,00001,00002... as de-multiplexing round indicator; \
   > 0000, 0001, 0002... as image number indicator; \
   > c1,c2,c3 or BF,CY5,CY3 or user_defined other names as image channel indicator; \
   > .tif_Results.xls as file extension.


### step 6, landmark cell and bead images, in ImageJ
 - INPUT: 
   > cell/ChannelGFP,TRITC,CY5 QUAD,BF_Seq0000_BF.tif \
   > bead/count00000_channelcy3,cy5,bf_seq0000_bf.tif 

 - RUN: 
   > measure the upper-left and bottom-right corner of the representative cell and bead image
 
 - OUTPUT:
   > UL_BR_{lane#}.xls (e.g. UL_BR_0.xls)
 
 - NOTE:
   > the rois in "UL_BR_{lane#}.xls" should be in the following order: cell.UL_well, cell.BR_well, bead.UL_well, bead.BR_well.
 

### step 7, bead optical demultiplexing, in Python
- INPUT: 

1. project_folder - the location of the project.
2. set_device 
   > n_lane: int - index of device lane \
   > total_lanes: int - number of lanes in the device. 1 for small, medium small device. 5 for large device. \
   > patches_per_lane: int - number of segments per device lane. 10 for medium small. 6 for large segment. \
   > first_border: int - x pixel position for the start of the first segment. \
   > border_gap: int - segment length in pixel along the x axis \
   > d_th: int - radius of well (pixel) - radius of bead (pixel). default 40; if use low resolution image for 
   demultiplexing, use 20.
3. set_demultiplexing
   > round_number: int - total number of demultiplexing cycles \
   > bead_channel: list - fluorescent channels used in demultiplexing [bf, S-cy5, Q-cy3] \
   > barcode_ref_fn: file name - bead optical barcode, sequence to binary code reference

- RUN:
```python
bead_intensity = BeadIntensity(project_folder)
# if use low resolution imaging setting
# bead_intensity.set_device(n_lane=0, total_lanes=1, patches_per_lane=10, first_border=200, border_gap=2470, d_th=15)
bead_intensity.set_device(n_lane=0, total_lanes=1, patches_per_lane=10, first_border=600, border_gap=7400, d_th=40)
bead_intensity.set_demultiplexing(bead_channel=['bf', 'cy5', 'cy3'], round_number=8, barcode_ref_fn = project_folder + 'bead/192_8mer_seq_reference.csv')
bead_intensity.generate_bead_intensity()
bead_intensity.obc_calling_correlation()
pickle_dump(bead_intensity, project_folder + 'bead_intensity.obj')
```

- OUTPUT: 

1. bead_intensity.obj
2. project/analysis 
   > probe_intensity.pdf \
   > obc_usage.pdf
 

### step 8, link sequencing data to imaging data, in Python

- INPUT:
1. set_bead_cell
   > bead_fn: BeadIntensity object file name, from step 7
   > bead_fn: CellIntensity object file name, from step 4
   > landmark_fn: file name - register markers of cell and bead scan

- RUN:
```python
cell_bead_link = CellBeadLink(project_folder)
cell_bead_link.set_bead_cell(bead_fn=project_folder + 'bead_intensity.obj', 
                             cell_fn=project_folder + 'cell_intensity.obj', 
                             landmark_fn=project_folder + 'well/UL_BR_0.xls')
cell_bead_link.link_obc_cell()
pickle_dump(cell_bead_link, project_folder + 'cell_bead_link.obj')
```

### step 9, link imaging to sequencing, in Python

- INPUT:
1. set_image
   > image_fn: file name - CellBeadLink object file
2. set_sequencing
   > cell_id: list, np.ndarray - cell index \
              file name - cell index \
   > use: category - 'index', use index column as cell index \
                     'columns', use columns names as cell index, 'gid', 'gene' are removed.
3. seq_image_link
   > replace: dict - {image patch: sequencing i7 index}

- RUN:
```python
image_seq_link = ImageSeqLink(project_folder)
image_seq_link.set_image(image_fn=project_folder + 'cell_bead_link.obj')
image_seq_link.set_sequencing(cell_id=project_folder + 'seq/PLZ003.edrops.matrix.hist.txt', use='index')
image_seq_link.seq_image_link(replace={0: 4, 1: 4, 2: 5, 3: 5, 4: 6, 5: 6, 6: 9, 7: 9, 8: 10, 9: 10})
pickle_dump(image_seq_link, project_folder + 'image_seq_link.obj')
```

- OUTPUT:
1. project/analysis
   > sc_image.matrix.linked.txt \
   > sc_seq.hist.linked.txt \
   > sc_seq.matrix.linked.txt

## objects

#### objects - bead intensity

The following code would generate a BeadIntensity object for each single lane, \
and store as an object using pickle. \
The BeadIntensity object contains: \
(1) self.bead_intensity_folder - the bead folder \
(2) self.channel_name - the channels of the image \
(3) self.n_lane - the number of lane, which the object related to \
(4) self.total_lanes - the total lanes of the device \
(5) self.n_patch_per_lane - the number of patches of each lane used for sequencing. \
(6) self.round - the number of rounds used for de-multiplexing \
(7) self.d_th - the distance threshold for register \
(8) self.bead_area - the area of the bead \
(9) self.bead_position - the location of the bead (pixel) \
(10) self.probe - the probe intensity of the beads \
(11) self.bg - the background intensity of the beads \
(12) self.probe_normalized - libear normalized probe intensity matrix \
(13) self.patch - the patch-id of the bead \
(14)self.obc_s - the s optical barcode of the bead \
(15)self.obc_q - the q optical barcode of the bead \
(16)self.round \
(17) self.obc - the optical barcode of the bead \
e.g. 0_19_23 means the 0th patch, 19th S-probe, 23th Q-probe \


#### objects - well_bead_cell link

The following code would generate a WellBeadCell object for each single lane, \
and store as an object using pickle. It uses well location as a reference, links the bead and the cell to their closest well. \
The WellBeadCell object contains: \
(1) self.n_lane - the number of lane, which the object related to \
(2) self.well_folder - the well folder \
(3) self.well_position - the well locations, as a reference for linking beads and cells \
(4) self.bead - the lane associated BeadIntensity object \
(5) self.rotation_matrix - the rotation matrix for aligning the bead and cell \
(6) self.bead_rotated_position - the rotated bead position \
(7) self.well_bead_link - bead id linked to the well \
(8) self.cell_folder - the cell folder \
(9) self.cell_position - the cell position \
(10) self.cell - the cell features, stored as a dict. key is the channel name, value is the measurement dataframe \
(11) self.well_cell_link - cell id linked to the well \
(12) self.obc_cell - the bead optical barcode linked to the cell id \
(13) self.sc_image_features_table - the fluorescent intensity measurement of all imaged single cells


#### objects - cell_features
The following code would generate a CellFeatures object for the whole device, \
and store as an object using pickle. \
The CellFeatures object contains: \
(1) self.cell_id - the cell ids \
(2) self.seqobc_cell - the sequencing cell id and the linked cell number id \
(3) self.cluster - the cluster of the cell \
(4) self.umap - the umap coordinate of the cell \
(5) self.count - the count matrix of the cell, cell id (row) by gene (col) \
(6) self.data - sequencing cell id and the well based fluorescent measurement \
(7) self.data_linked - sequencing cell id and the well based fluorescent measurement, with only linked single cells \
(8) self.sc_image_matrix - image feature matrix. sequencing cell id and the single cell based fluorescent measurement, with only linked single cells \
(9) self.sc_seq_matrix - count matrix. sequencing cell id and mRNA count with only linked single cells \
(10) self.sc_seq_hist - count hist matrix. sequencing cell id and mRNA count summary with only linked single cells

