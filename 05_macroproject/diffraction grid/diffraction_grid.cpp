#include <dolfin.h>
#include <string>
#include "Wave_equation.h"
#include "Square.h"
#include "Add.h"
#include "Mean.h"
#include <iostream>

using namespace dolfin;

class Boundary_Func : public Expression // Функция на границе
{
public:
  Boundary_Func() : t(0) {}
  // Define boundary condition
  void eval(Array<double> &values, const Array<double> &x) const
  {
      values[0] = 30 * cos(3000 * t);
  }
  double t;
};

class DirichletBoundary : public SubDomain // граница
{
  bool wavering_zones(const Array<double> &x, double b, int N) const
  {
    bool indicator = false;
    for (int i = 0; i < N; i++)
      if ((x[1] > ((-(2 * N + 1) / 2 + 2 * i) * b)) && (x[1] < ((-(2 * N + 1) / 2 + (2 * i + 1)) * b)))
        indicator = true;
    return indicator;
  }

  bool inside(const Array<double> &x, bool on_boundary) const
  {
    double c = 0.005;
    return ((x[0] < c) && wavering_zones(x, 0.02, 25));
  }
};

class Indicator : public Expression // экран
{
  void eval(Array<double> &values, const Array<double> &x) const
  {
    double E = 0.03; // трение
    double c = 0.03;
    // double top = 1; //верхняя грань по y
    // double b = 0.05; //размер щели; d - период, d = 2*b
    // if (((x[1] > 9/2*b) or (x[1] < -9/2*b) or ((x[1] < 7*b/2) and (x[1] > 5*b/2)) or ((x[1] < 3*b/2) and (x[1] > b/2)) or
    // ((x[1] < -b/2) and (x[1] > -3*b/2)) or ((x[1] < -5*b/2) and (x[1] > -7*b/2))) and ((x[0] > -0.7) and (x[0] < -0.65))) values[0] = c;
    // else if ((x[0] > -0.7) && ((x[1] > top-0.1) || (x[1]) < -top+0.1)) values[0] = c;
    // else if (x[0] > 0.85) values[0] = c;
    // //else if (x[0] < 0) values[0] = 2*x[0] + 1;
    // else values[0] = 0;
    // if (x[0] > 1.5-c) values[0] = E;
    // else if ((x[1] > x[0] + 1 - c) || (x[1] < -x[0] - 1 + c)) values[0] = E;
    // else values[0] = 0;
    values[0]= 0;
  }
};

int main()
{
  // std::array<Point, 2> a1 = {Point(-1, -1), Point(1, 1)};
  // std::array<long unsigned int, 2> a2 = {350, 500};
  auto mesh = std::make_shared<Mesh>("mesh3.xml");
  auto V = std::make_shared<Wave_equation::FunctionSpace>(mesh);

  auto u0 = std::make_shared<Boundary_Func>();
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc(V, u0, boundary);

  double dt = 0.00005;

  auto u_pr = std::make_shared<Function>(V);
  auto u_prpr = std::make_shared<Function>(V);
  auto k = std::make_shared<Constant>(dt);
  auto ind = std::make_shared<Indicator>();

  Wave_equation::BilinearForm a(V, V);
  Wave_equation::LinearForm L(V);
  Square::BilinearForm a_sq(V, V);
  Square::LinearForm L_sq(V);
  Add::BilinearForm a_add(V, V);
  Add::LinearForm L_add(V);

  auto u = std::make_shared<Function>(V);
  auto I = std::make_shared<Function>(V);     // интенсивность, вычисленная на этом шаге
  auto I_pr = std::make_shared<Function>(V);  // интенсивность на предыдущем шаге
  auto I_all = std::make_shared<Function>(V); // интенсивность сейчас
  a.dt = k;
  L.dt = k;
  a.indicator = ind;
  L.indicator = ind;

  int steps = 100000;
  for (int i = 0; i < steps; i++)
  {
    u0->t = i * dt;
    DirichletBC bc(V, u0, boundary);
    L.u_pr = u_pr;
    L.u_prpr = u_prpr;
    //  L.g = g;
    solve(a == L, *u, bc); //основное вычисление

    //L_sq.u = u;
    //solve(a_sq == L_sq, *I);
    std::cout << u0->t << '\n';

    L_add.I_pr = I_pr;
    L_add.E = u;
    solve(a_add == L_add, *I_all);
    *I_pr = *I_all;
    if (i % 25 == 0)
    {
      std::string fileName1 = "snapshots_diffraction_grid_N=25_omega=3000/wave_equation-" + std::to_string(i / 25) + ".pvd";
      File file1(fileName1);
      file1 << *u;
    }
    if (i % 125 == 0)
    {
      std::string fileName = "snapshots_diffraction_grid_I_N=25_omega=3000/wave_equation-" + std::to_string(i / 125) + ".pvd"; // строит суммы квадратов напряжённостей поля
      File file(fileName);
      file << *I_all;
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
  std::string fileName2 = "snapshots_Yong's_experiment/I.pvd"; // Строит именно распределение интенсивности
  File file2(fileName2);
  file2 << *I_final;

  return 0;
}