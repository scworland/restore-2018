
# set working directory
import os,sys
os.chdir(os.path.expanduser('~/Documents/Restore'))

def dbf2feather():
   
    import glob
    path = "data/basinchars/geodatabases/*.dbf"
    for fname in glob.glob(path):

        import pandas as pd
        from osgeo import ogr
        import re

        t = ogr.Open(fname)
        lyr = t.GetLayer(0)
        cols = {}

        lyrdef = lyr.GetLayerDefn()
        for i in range(lyrdef.GetFieldCount()):
            fd = lyrdef.GetFieldDefn(i)
            fn = fd.name
            cols[fn] = []

        for feat in lyr:
            for k in cols.keys():
                cols[k].append(feat.GetField(k))

        df = pd.DataFrame.from_dict(cols)
        fname_feather=fname.replace('data/basinchars/geodatabases/','').replace('.dbf','.feather')

        df.to_feather(fname_feather)



dbf2feather()