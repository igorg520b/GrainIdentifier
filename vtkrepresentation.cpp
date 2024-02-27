#include "vtkrepresentation.h"

VtkRepresentation::VtkRepresentation()
{
    int nLut = sizeof lutArrayPastel / sizeof lutArrayPastel[0];
    hueLut_pastel->SetNumberOfTableValues(nLut);
    for ( int i=0; i<nLut; i++)
        hueLut_pastel->SetTableValue(i, lutArrayPastel[i][0],
                                     lutArrayPastel[i][1],
                                     lutArrayPastel[i][2], 1.0);
    hueLut_pastel->SetTableRange(0,39);

    visualized_values->SetName("visualized_values");

    ugrid->SetPoints(points_mesh);
    dataSetMapper->SetInputData(ugrid);
    actor_grains->SetMapper(dataSetMapper);
    actor_grains->GetProperty()->SetOpacity(0.5);


    ugrid->GetCellData()->AddArray(visualized_values);
    ugrid->GetCellData()->SetActiveScalars("visualized_values");
    dataSetMapper->SetLookupTable(hueLut_pastel);
    dataSetMapper->UseLookupTableScalarRangeOn();
    dataSetMapper->SetScalarModeToUseCellData();
    dataSetMapper->ScalarVisibilityOn();

    source_cube->SetCenter(0.5,0.5,0.5);
    mapper_cube->SetInputConnection(source_cube->GetOutputPort());
    actor_cube->SetMapper(mapper_cube);
    actor_cube->GetProperty()->SetRepresentationToWireframe();
    actor_cube->GetProperty()->SetColor(0.1,0.1,0.1);
}


void VtkRepresentation::SynchronizeTopology()
{
    points_mesh->SetNumberOfPoints(gp->vertices2.size());

    for(int i=0;i<gp->vertices2.size();i++)
    {
        Eigen::Vector3f &v = gp->vertices2[i];
        double x[3] {v[0],v[1],v[2]};
        points_mesh->SetPoint((vtkIdType)i, x);
    }

    cellArray->Reset();


    visualized_values->SetNumberOfValues(gp->elems2.size());

    for(int i=0;i<gp->elems2.size();i++)
    {
        vtkIdType pts[4];
        for(int k=0;k<4;k++) pts[k] = gp->elems2[i][k];
        cellArray->InsertNextCell(4, pts);
        visualized_values->SetValue(i, gp->elems2[i][4]%40);
    }
    ugrid->SetCells(VTK_TETRA, cellArray);


}


