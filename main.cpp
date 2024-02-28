#include "mainwindow.h"

#include <QApplication>
#include <gmsh.h>
#include <omp.h>
#include <spdlog/spdlog.h>

int main(int argc, char *argv[])
{
    gmsh::initialize();
    spdlog::info("testing threads {}", omp_get_max_threads());
#pragma omp parallel
    {     spdlog::info("{}", omp_get_thread_num()); }
    std::cout << std::endl;

    QApplication a(argc, argv);
    MainWindow w;
//    w.show();
    w.showMaximized();
    return a.exec();
}
