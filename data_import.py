import pandas as pd


def extract_as_object(name):
    columns=['y-span', 'Chord', 'Ai', 'Cl', 'PCd', 'ICd', 'CmGeom', 'CmAirf@chord/4', 'XTrtop', 'XTrBot', 'XCP', 'BM']
    class span_pos:
        def __init__(self, span, chord, Ai, Cl, lCd, Cm):
            self.span = span # span-wise location
            self.chord = chord
            self.Ai = Ai # AOA, induced
            self.Cl = Cl # Cl
            self.lCd = lCd # induced cd
            self.Cm = Cm # quarter chord moment

    data=pd.read_csv(name,sep=r'\s+',header=None,names=columns)
    data = data.drop(0).reset_index(drop=True)
    useful_data=data[['y-span','Chord','Ai','Cl','ICd','CmGeom']]
    num_rows = len(useful_data)
    col_count=useful_data.shape[1]
    data_points=[]
    for j in range(num_rows):
        singular_point=[]
        for i in range(col_count):
            singular_point.append(float(useful_data.iloc[j, i]))
        data_points.append(singular_point)
    span_objects=[]
    for i in data_points:
        span_obj=span_pos(i[0], i[1], i[2], i[3], i[4], i[5])
        span_objects.append(span_obj)
    return span_objects

def extract_as_list(name):
    columns=['y-span', 'Chord', 'Ai', 'Cl', 'PCd', 'ICd', 'CmGeom', 'CmAirf@chord/4', 'XTrtop', 'XTrBot', 'XCP', 'BM']
    data=pd.read_csv(name,sep=r'\s+',header=None,names=columns)
    data = data.drop(0).reset_index(drop=True)
    useful_data=data[['y-span','Chord','Ai','Cl','ICd','CmGeom']]
    num_rows = len(useful_data)
    col_count=useful_data.shape[1]
    data_points=[]
    for j in range(col_count):
        singular_point=[]
        for i in range(num_rows):
            singular_point.append(float(useful_data.iloc[i, j]))
        data_points.append(singular_point)
    return data_points

