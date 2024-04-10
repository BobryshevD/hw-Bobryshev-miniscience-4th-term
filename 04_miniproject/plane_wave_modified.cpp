#include <dolfin.h>
#include <string>
#include "Wave_equation.h"
#include "Square.h"
#include "Add.h"
#include "Mean.h"
#include <iostream>

using namespace dolfin;

class Boundary_Func : public Expression
{
public:
  Boundary_Func() : t(0) {}
  // Define boundary condition
  void eval(Array<double> &values, const Array<double> &x) const
  {
    double c = 1.0 / 40;
    if (x[0] < c - 1)
      values[0] = 30 * cos(100 * t);
    else
      values[0] = 0;
  }
  double t;
};

class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double> &x, bool on_boundary) const
  {
    double c = 1.0 / 50;
    return (x[0] < c - 1) || ((x[0] > -c-0.5) && (x[0] < c-0.5) && ((x[1] > 5*c) || (x[1] < -5*c)));
  }
};

int main()
{
  std::array<Point, 2> a1 = {Point(-1, -1), Point(3, 1)};
  std::array<long unsigned int, 2> a2 = {600, 400};
  auto mesh = std::make_shared<Mesh>(RectangleMesh::create(a1, a2, CellType::Type::triangle));
  auto V = std::make_shared<Wave_equation::FunctionSpace>(mesh);

  auto u0 = std::make_shared<Boundary_Func>();
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc(V, u0, boundary);

  double dt = 0.0005;

  auto u_pr = std::make_shared<Function>(V);
  auto u_prpr = std::make_shared<Function>(V);
  auto k = std::make_shared<Constant>(dt);

  Wave_equation::BilinearForm a(V, V);
  Wave_equation::LinearForm L(V);
  Square::BilinearForm a_sq(V, V);
  Square::LinearForm L_sq(V);
  Add::BilinearForm a_add(V, V);
  Add::LinearForm L_add(V);

  auto u = std::make_shared<Function>(V);
  auto g = std::make_shared<Constant>(0);
  auto I = std::make_shared<Function>(V);     // интенсивность, вычисленная на этом шаге
  auto I_pr = std::make_shared<Function>(V);  // интенсивность на предыдущем шаге
  auto I_all = std::make_shared<Function>(V); // интенсивность сейчас
  a.dt = k;
  int steps = 3000;
  for (int i = 0; i < steps; i++)
  {
    u0->t = i * dt;
    DirichletBC bc(V, u0, boundary);
    L.u_pr = u_pr;
    L.u_prpr = u_prpr;
  //  L.g = g;
    solve(a == L, *u, bc);
    L_sq.u = u;
    solve(a_sq == L_sq, *I);
    std::cout << u0->t << '\n';

    L_add.I_pr = I_pr;
    L_add.I_all = I;
    solve(a_add == L_add, *I_all);
    *I_pr = *I_all;
    if (i%5 == 0)
    {
      std::string fileName = "snapshots_plane_wave_modified_I/wave_equation-" + std::to_string(i/5) + ".pvd";
      File file(fileName);
      file << *I_all;
      std::string fileName1 = "snapshots_plane_wave_modified/wave_equation-" + std::to_string(i/5) + ".pvd";
      File file1(fileName1);
      file1 << *u;
    }

    *u_prpr = *u_pr;
    *u_pr = *u;
  }

  auto I_final = std::make_shared<Function>(V);
  Mean::BilinearForm a_mean(V, V);
  Mean::LinearForm L_mean(V);
  auto i = std::make_shared<Constant>(steps);
  L_mean.i = i;
  L_mean.I_all = I_all;
  solve(a_mean == L_mean, *I_final);
  std::string fileName2 = "snapshots_plane_wave_modified/I.pvd";
  File file2(fileName2);
  file2 << *I_final;

  return 0;
}