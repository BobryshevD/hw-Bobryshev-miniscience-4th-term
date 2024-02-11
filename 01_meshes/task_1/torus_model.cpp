#include <set>
#include <gmsh.h>


int main(int argc, char **argv){
    gmsh::initialize();

    gmsh::model::add("t2");

    double lc = 1e-1;
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
    gmsh::model::geo::addPoint(0, 2, 0, lc, 2);
    gmsh::model::geo::addPoint(0, 4, 0, lc, 3);//для большой окружности в плоскости xy

    gmsh::model::geo::addPoint(0, 1, 0, lc, 4);//для малой окружности в плоскости xy
    gmsh::model::geo::addPoint(0, 3, 0, lc, 5);

    gmsh::model::geo::addPoint(0, 0.5, 0, lc, 10);
    gmsh::model::geo::addPoint(0, 3.5, 0, lc, 11);

    gmsh::model::geo::addPoint(0, 0.25, 0, lc, 12);//для большой окружности в плоскости xy
    gmsh::model::geo::addPoint(0, 3.75, 0, lc, 13);

    gmsh::model::geo::addPoint(0, 0.75, 0, lc, 14);//для малой окружности в плоскости xy
    gmsh::model::geo::addPoint(0, 3.25, 0, lc, 15);


    gmsh::model::geo::addCircleArc(1, 2, 3, 1);//первая окружность в xy
    gmsh::model::geo::addCircleArc(3, 2, 1, 2);

    gmsh::model::geo::addCircleArc(4, 2, 5, 3);//вторая окружность в xy
    gmsh::model::geo::addCircleArc(5, 2, 4, 4);

    gmsh::model::geo::addCircleArc(4, 10, 1, 5, 1, 0, 0);//первая окружность в xz
    gmsh::model::geo::addCircleArc(1, 10, 4, 6, 1, 0, 0);

    gmsh::model::geo::addCircleArc(3, 11, 5, 7, 1, 0, 0);//вторая окружность в xz
    gmsh::model::geo::addCircleArc(5, 11, 3, 8, 1, 0, 0);


    gmsh::model::geo::addCircleArc(12, 2, 13, 9);
    gmsh::model::geo::addCircleArc(13, 2, 12, 10);

    gmsh::model::geo::addCircleArc(14, 2, 15, 11);
    gmsh::model::geo::addCircleArc(15, 2, 14, 12);

    gmsh::model::geo::addCircleArc(13, 11, 15, 13, 1, 0, 0);
    gmsh::model::geo::addCircleArc(15, 11, 13, 14, 1, 0, 0);

    gmsh::model::geo::addCircleArc(14, 10, 12, 15, 1, 0, 0);
    gmsh::model::geo::addCircleArc(12, 10, 14, 16, 1, 0, 0); //просто нарисовали два вложенных тора


    gmsh::model::geo::addCurveLoop({1, 7, -3, 5}, 1);
    gmsh::model::geo::addSurfaceFilling({1}, 1);

    gmsh::model::geo::addCurveLoop({8, -1, 6, 3}, 2);
    gmsh::model::geo::addSurfaceFilling({2}, 2);

    gmsh::model::geo::addCurveLoop({7, 4, 5, -2}, 3);
    gmsh::model::geo::addSurfaceFilling({3}, 3);

    gmsh::model::geo::addCurveLoop({8, 2, 6, -4}, 4);
    gmsh::model::geo::addSurfaceFilling({4}, 4);

    gmsh::model::geo::addCurveLoop({15, -10, 13, 12}, 5);
    gmsh::model::geo::addSurfaceFilling({5}, 5);

    gmsh::model::geo::addCurveLoop({16, -12, 14, 10}, 6);
    gmsh::model::geo::addSurfaceFilling({6}, 6);

    gmsh::model::geo::addCurveLoop({15, 9, 13, -11}, 7);
    gmsh::model::geo::addSurfaceFilling({7}, 7);
    
    gmsh::model::geo::addCurveLoop({16, 11, 14, -9}, 8);
    gmsh::model::geo::addSurfaceFilling({8}, 8);    


    gmsh::model::geo::addSurfaceLoop({1, 2, 3, 4}, 1); 
    gmsh::model::geo::addSurfaceLoop({5, 6, 7, 8}, 2);

     gmsh::model::geo::addVolume({1, 2});

    gmsh::model::geo::synchronize();

    gmsh::model::mesh::generate(3);


    gmsh::write("torus_model.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize(); 
}  

