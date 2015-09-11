#!/usr/bin/env python

#
# Class definition for object to hold and query variation data.
# Kristian Jensen, Januar 2015
#


class DataTable(list):
    '''A class for holding records of data. Provides method for filtering and summarizing the data.
    Initialize with an iterable of sequences and optionally a sequence of column names.
    '''
    def __init__(self,data,columns=None):
        if columns is None:
            columns = range(len(data[0]))
        assert len(columns) == len(set(columns)) # Assert column uniqueness
        self.columns = tuple(_ for _ in columns)
        for record in data:
            assert len(record) == len(columns) # Assert correct and equal length of records
            self.append(record)
    
    def columnSet(self,column_name):
        "Returns a set containing the values in the given column"
        index = self.columns.index(column_name)
        name_list = []
        for record in self:
            name_list.append(record[index])
        return set(name_list)
    
    def subset(self,column,value):
        "Returns a new DataTable containing only the records where <column> == <value>"
        index = self.columns.index(column)
        new_data = []
        for record in self:
            if record[index] == value:
                new_data.append(record)
        return DataTable(new_data,self.columns)
        

if __name__ == "__main__":
    
    raw_data = open("/Users/kristianjensen/Documents/DTU/speciale/Data/Uniprot_natural_variations.txt","r").read().splitlines()
    data = [_.split("\t") for _ in raw_data]
    columns = ["gene","start","end","strain","original","variant"]
    
    #test = [["hej","med"],["dig","der"],["du","der"]]
    #columns = ["first","second"]
    
    obj = DataTable(data,columns)
    
    print obj.columns
    
    for strain in sorted(obj.columnSet("strain")):
        print strain 
        