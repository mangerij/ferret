from paraview import *


reader = GetActiveSource()
# SVReader(DetectNumericColumns = True, FieldDelimiterCharacters = ",", HaveHeaders = True, FileName = "/home/qkhan/work/scalfmm/Build/20k10z.4.csv")

#classRef = CSVReader()

#if not isinstance(reader, classRef.__class__):
#    exit(-1)

filename = reader.FileName[0]
print filename
nbZones = int(filename.split('.')[0].split('_')[-1][0:-1])

selection = SelectionQuerySource(FieldType = "ROW", QueryString = "zone >= 0")
extractor = ExtractSelection(Input = reader, Selection = selection)
points = TableToPoints(Input = extractor, XColumn = "x", YColumn = "y", ZColumn = "z")

repr = GetRepresentation()
repr.ColorArrayName = 'zone'
repr.LookupTable = AssignLookupTable(points.PointData['zone'], "Cool to Warm")
Show()

#for i in range(nbZones):

