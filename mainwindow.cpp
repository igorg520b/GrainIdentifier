#include "mainwindow.h"
#include "./ui_mainwindow.h"

#include <filesystem>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);


    representation.gp = &gp;

    // VTK
    qt_vtk_widget = new QVTKOpenGLNativeWidget();
    qt_vtk_widget->setRenderWindow(renderWindow);
    renderer->SetBackground(1.0,1.0,.9);
    renderWindow->AddRenderer(renderer);
    setCentralWidget(qt_vtk_widget);

//    renderer->AddActor(representation.actor_grains);
    renderer->AddActor(representation.actor_cube);
    renderer->AddActor(representation.actor_MPM);

    std::string fileName = "/home/s2/Documents/neper_data/1k.msh";
    gp.LoadMSH(fileName);

    std::string pointsFile = "/home/s2/Projects-other/2m.h5";
    gp.Load_Points_HDF5(pointsFile);

    gp.IdentifyGrains();
    representation.SynchronizeTopology();
}

MainWindow::~MainWindow()
{
    delete ui;
}
