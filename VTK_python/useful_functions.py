import vtk
import numpy as np
import math
import collections
from vtk.util.numpy_support import vtk_to_numpy
def read_vtu(filename):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def read_vtp(filename):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()
def read_read(filename):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader

def get_normals(poly_data):
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(poly_data.GetOutputPort())
    normals.ComputePointNormalsOn()
    normals.ComputeCellNormalsOff()
    normals.Update()
    pd = normals.GetOutput()
    pd.GetPointData().SetActiveScalars('Normals')
    norms = vtk_to_numpy(pd.GetPointData().GetArray(1))
    return norms

def get_curvature_mean (poly_data):
    curvatures = vtk.vtkCurvatures()
    curvatures.SetInputConnection(poly_data.GetOutputPort())
    curvatures.SetCurvatureTypeToMean()
    curvatures.Update()
    pd = curvatures.GetOutput()
    pd.GetPointData().SetActiveScalars('Mean_Curvature')
    Mean = vtk_to_numpy(pd.GetPointData().GetArray(1))
    return Mean

def get_curvature_gaussian (poly_data):
    curv = vtk.vtkCurvatures()
    curv.SetInputConnection(poly_data.GetOutputPort())
    curv.SetCurvatureTypeToGaussian()
    curv.Update()
    pd = curv.GetOutput()
    pd.GetPointData().SetActiveScalars('Gauss_Curvature')
    gaussian = vtk_to_numpy(pd.GetPointData().GetArray(1))
    return gaussian

def get_neighbors_around_minima (f,polyData,rd):
    points = np.array(polyData.GetPoints().GetData())
    gaussian = get_curvature_gaussian (rd)
    Mean = get_curvature_mean(rd)
    normals = get_normals(rd)
    T_Tree = vtk.vtkKdTreePointLocator()
    T_Tree.SetDataSet(polyData)
    T_Tree.BuildLocator()
    neighbours = []
    neighbours.append([points[f][0],points[f][1],points[f][2],0,gaussian[f],Mean[f],normals[f][0],normals[f][1],normals[f][2]])
    f = [points[f][0],points[f][1],points[f][2]]
    result = vtk.vtkIdList()
    num_of_neighbors = 200
    T_Tree.FindClosestNPoints(num_of_neighbors, f, result);
    for i in range (0,result.GetNumberOfIds()):
        point_ind = result.GetId(i)
        closestPoint = [0.0, 0.0, 0.0]
        T_Tree.GetDataSet().GetPoint(point_ind, closestPoint)
        distance = float(math.sqrt(vtk.vtkMath.Distance2BetweenPoints(closestPoint, f)))
        neighbours.append([closestPoint[0], closestPoint[1], closestPoint[2],distance,gaussian[point_ind],Mean[point_ind],normals[point_ind][0],normals[point_ind][1],normals[point_ind][2]])
        y=np.array([np.array(xi) for xi in neighbours])
    return np.array(neighbours)       
def visualize_data(model):
    scalar_range = model.GetScalarRange()

    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(model)
    mapper.SetScalarRange(scalar_range)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(1, 1, 1)  
    renderer_window = vtk.vtkRenderWindow()
    renderer_window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderer_window)
    interactor.Initialize()
    interactor.Start()
    return renderer_window

def write_vtp(file_name, polyData):
    out = vtk.vtkXMLPolyDataWriter()
    out.SetInputDataObject(polyData)
    out.SetDataModeToAscii()
    out.SetFileName(file_name)
    out.Write()
    return 1

def write_vtu_ascii(file_name, polyData):
    out = vtk.vtkXMLUnstructuredGridWriter()
    out.SetInputDataObject(polyData)
    out.SetDataModeToAscii()
    out.SetFileName(file_name)
    out.Write()
    return     


def get_minimum(polyData):
    num_of_points = polyData.GetNumberOfPoints()
    #print(num_of_points)
    points = np.array(polyData.GetPoints().GetData())
    #used getArray(0) in case there are more than one array
    Safty_factor = np.array(polyData.GetPointData().GetArray(0))
    point_cells = vtk.vtkIdList();
    cell_points = vtk.vtkIdList();
    x=0
    local_minima=[]
    Safe_fac=[]
    for point in range (0 , num_of_points):
        minima = 1
        polyData.GetPointCells( point, point_cells )
        num_of_cells = point_cells.GetNumberOfIds();
        for cell in range (0 , num_of_cells):
            neighbor_cell = point_cells.GetId( cell )
            polyData.GetCellPoints(neighbor_cell , cell_points )
            num_cell_points = cell_points.GetNumberOfIds();
            for neighbor_point in range (0, num_cell_points) :
                neighbor_point_id = cell_points.GetId(neighbor_point)
                if ( Safty_factor[neighbor_point_id] <= Safty_factor[point] and neighbor_point_id != point):
                    minima = 0
        if ( minima == 1 and Safty_factor[point] < 40):
            local_minima.append(point)
            #p_location = polyData.GetPoint( point )
            #print(Safty_factor[point])
            #print(p_location[0] ,p_location[1], p_location[2])
    dic_ind_Sf=dict(zip(local_minima,Safty_factor[local_minima]))
    sorted_ = collections.OrderedDict(sorted(dic_ind_Sf.items(), key=lambda x: x[1]))
    return sorted_

