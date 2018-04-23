#!/usr/bin/env python

import pandas as pd
import numpy as np
coord = pd.read_table('mOTU.v2b.padded.coords_for_centroids', header = None) #mOTU.v2b.centroids.reformatted.padded.coord

index = coord.index
columns = ['gene_id','external_id','sequence_id',
		   'type','gene_info','length',
		   'start','end','strand','start_codon',
		   'stop_codon','gc']

annotations = pd.DataFrame(index=index, columns=columns)
annotations.loc[:,'gene_id'] = coord.index + 1
annotations.loc[:,'external_id'] = coord.iloc[:,0]
annotations.loc[:,'sequence_id'] = coord.iloc[:,0]
annotations.loc[:,'type'] = 'CDS'
annotations.loc[:,'gene_info'] = '<annotation product= NaN />'
annotations.loc[:,'length'] = coord.iloc[:,3] - coord.iloc[:,2] +1
annotations.loc[:,'start'] = coord.iloc[:,2]
annotations.loc[:,'end'] = coord.iloc[:,3]
annotations.loc[:,'strand'] = '+'
annotations.loc[:,'start codon'] = '0'
annotations.loc[:,'end codon'] = '0'
annotations.loc[:,'GC content'] = 'NaN'


annotations.to_csv('mOTU.v2b.centroids.reformatted.padded.annotations', index = False, header = False, sep ='\t')
