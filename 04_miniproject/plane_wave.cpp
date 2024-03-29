#include <dolfin.h>
#include <string>
#include "Wave_equation.h"
#include <iostream>

using namespace dolfin;

class Boundary_Func : public Expression
{
public:
  Boundary_Func(): t(0) {}
  // Define boundary condition
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 30*cos(50*t);
  }
  double t;
};

class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    double c = 1.0/40;
    return (x[0] < c-1);
  }
};

int main()
{
  std::array<Point, 2> a1 = {Point(-1, -1), Point(1, 1)};
  std::array<long unsigned int, 2> a2 = {80, 80};
  auto mesh = std::make_shared<Mesh>(RectangleMesh::create(a1, a2, CellType::Type::triangle));
  auto V = std::make_shared<Wave_equation::FunctionSpace>(mesh);
  
  auto u0 = std::make_shared<Boundary_Func>();
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc(V, u0, boundary);

  double dt = 0.001;

  auto u_pr = std::make_shared<Function>(V);
  auto u_prpr = std::make_shared<Function>(V);
  auto k = std::make_shared<Constant>(dt);

  Wave_equation::BilinearForm a(V, V);
  Wave_equation::LinearForm L(V);
  auto u = std::make_shared<Function>(V);
  auto g = std::make_shared<Constant>(0);
  a.dt = k;

  for (int i = 0; i<1000; i++)
  {
    u0->t = i*dt;
    DirichletBC bc(V, u0, boundary);
    L.u_pr = u_pr;
    L.u_prpr = u_prpr;
    solve(a == L, *u, bc);
    
    std::cout << u0->t << '\n';

    std::string fileName = "snapshots_plane_wave/wave_equation-" + std::to_string(i) + ".pvd";
    File file(fileName);
    file << *u;

    *u_prpr = *u_pr;
    *u_pr = *u;
  }


  return 0;
}
