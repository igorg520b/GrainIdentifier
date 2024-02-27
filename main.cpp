#include "mainwindow.h"

#include <QApplication>
#include <gmsh.h>

int main(int argc, char *argv[])
{
    gmsh::initialize();

    QApplication a(argc, argv);
    MainWindow w;
//    w.show();
    w.showMaximized();
    return a.exec();
}
