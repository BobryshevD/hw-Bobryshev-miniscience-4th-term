#include <set>
#include <cmath>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("t13");

  try {
    gmsh::merge("d20.stl");
  } catch(...) {
    gmsh::logger::write("Could not load STL mesh: bye!");
    gmsh::finalize();
    return 0;
  }

  double angle = 40;

//  bool forceParametrizablePatches = false;
  bool includeBoundary = true;
  double curveAngle = 180;

  gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                      curveAngle * M_PI / 180.);


  gmsh::model::mesh::createGeometry();
  std::vector<std::pair<int, int> > s;
  gmsh::model::getEntities(s, 2);
  std::vector<int> sl;
  for(auto surf : s) sl.push_back(surf.second);
  int l = gmsh::model::geo::addSurfaceLoop(sl);
  gmsh::model::geo::addVolume({l});

  gmsh::model::geo::synchronize();

  // We specify element sizes imposed by a size field, just because we can :-)
  bool funny = true; // true;
  int f = gmsh::model::mesh::field::add("MathEval");
  if(funny)
    gmsh::model::mesh::field::setString(f, "F", "2*Sin((x+y)/5) + 3");
  else
    gmsh::model::mesh::field::setString(f, "F", "4");
  gmsh::model::mesh::field::setAsBackgroundMesh(f);

  gmsh::model::mesh::generate(2);

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();
  return 0;
}