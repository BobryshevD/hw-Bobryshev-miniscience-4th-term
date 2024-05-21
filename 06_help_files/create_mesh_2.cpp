#include <set>
#include <gmsh.h>


int main(int argc, char **argv){
    gmsh::initialize();

    gmsh::model::add("t2");

    double lc = 0.008;
    double lc2 = 0.1;
    //gmsh::model::geo::addPoint(0, -1, 0, lc, 1);
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
    //gmsh::model::geo::addPoint(0, 1, 0, lc, 3);
    gmsh::model::geo::addPoint(0, 1, 0, lc, 2);
    gmsh::model::geo::addPoint(1, 0, 0, lc, 3);
    gmsh::model::geo::addPoint(1, 1, 0, lc, 4);
    gmsh::model::geo::addPoint(3, 0, 0, lc2, 5);
    gmsh::model::geo::addPoint(3, 1, 0, lc2, 6);
    
    //gmsh::model::geo::addLine(1, 2);
    gmsh::model::geo::addLine(1, 3, 1);
    gmsh::model::geo::addLine(3, 4, 2);
    gmsh::model::geo::addLine(4, 2, 3);
    //gmsh::model::geo::addLine(4, 3);
    //gmsh::model::geo::addLine(3, 1);
    gmsh::model::geo::addLine(2, 1, 4);

    gmsh::model::geo::addLine(4, 6, 5);
    gmsh::model::geo::addLine(6, 5, 6);
    gmsh::model::geo::addLine(5, 3, 7);
    

    gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);    
    gmsh::model::geo::addCurveLoop({2, 5, 6, 7}, 2);
    gmsh::model::geo::addPlaneSurface({2}, 2);

    gmsh::model::geo::synchronize();

    gmsh::model::mesh::generate(2);


    gmsh::write("test_mesh.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize(); 
}  

