#include <set>
#include <gmsh.h>


int main(int argc, char **argv){
    gmsh::initialize();

    gmsh::model::add("t2");

    double lc = 0.0035;
    double lc2 = 0.5;
    gmsh::model::geo::addPoint(0, -1, 0, lc, 1);
    gmsh::model::geo::addPoint(0.5, -1, 0, lc, 2);
    gmsh::model::geo::addPoint(0, 1, 0, lc, 3);
    gmsh::model::geo::addPoint(0.5, 1, 0, lc, 4);
    gmsh::model::geo::addPoint(1.7, 2.2, 0, lc, 5);
    gmsh::model::geo::addPoint(1.7, -2.2, 0, lc, 6);
    gmsh::model::geo::addPoint(0, 5, 0, lc2, 7);
    gmsh::model::geo::addPoint(1.7, 5, 0, lc2, 8);
    gmsh::model::geo::addPoint(7, 5, 0, lc2, 9);
    gmsh::model::geo::addPoint(0, -5, 0, lc2, 10);
    gmsh::model::geo::addPoint(1.7, -5, 0, lc2, 11);
    gmsh::model::geo::addPoint(7, -5, 0, lc2, 12);
    


    gmsh::model::geo::addLine(1, 2, 1);
    gmsh::model::geo::addLine(2, 6, 2);
    gmsh::model::geo::addLine(6, 5, 3);
    gmsh::model::geo::addLine(5, 4, 4);
    gmsh::model::geo::addLine(4, 3, 5);
    gmsh::model::geo::addLine(3, 1, 6);

    gmsh::model::geo::addLine(3, 7, 7);
    gmsh::model::geo::addLine(7, 8, 8);
    gmsh::model::geo::addLine(8, 9, 9);
    gmsh::model::geo::addLine(9, 12, 10);
    gmsh::model::geo::addLine(12, 11, 11);
    gmsh::model::geo::addLine(11, 10, 12);
    gmsh::model::geo::addLine(10, 1, 13);
    gmsh::model::geo::addLine(1, 2, 14);
    gmsh::model::geo::addLine(2, 6, 15);
    gmsh::model::geo::addLine(6, 5, 16);
    gmsh::model::geo::addLine(5, 4, 17);
    gmsh::model::geo::addLine(4, 3, 18);
    


    gmsh::model::geo::addCurveLoop({1, 2, 3, 4, 5, 6}, 1);
    gmsh::model::geo::addCurveLoop({7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}, 2);

    gmsh::model::geo::addPlaneSurface({1}, 1);    
    gmsh::model::geo::addPlaneSurface({2}, 2);

    gmsh::model::geo::synchronize();

    gmsh::model::mesh::generate(2);


    gmsh::write("mesh.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize(); 
}  

