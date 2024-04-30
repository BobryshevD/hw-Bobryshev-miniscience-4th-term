#include <dolfin.h>
#include <string>
#include "Wave_equation.h"
#include <iostream>

using namespace dolfin;

class Boundary_Func : public Expression  //функция на границе
{
public:
    Boundary_Func() : t(0) {}
    // Define boundary condition
    void eval(Array<double> &values, const Array<double> &x) const
    {
        if (t < 0.04) values[0] = 30 * sin(300 * t);
        else values[0] = 0;   
    }
    double t;
};

class DirichletBoundary : public SubDomain //граница
{
    bool inside(const Array<double> &x, bool on_boundary) const
    {
        double c = 1.0 / 40;
        return (x[0] < c - 1);
    }
};

class Indicator : public Expression //!!!Область с другим показателем преломления!!!
{
    void eval(Array<double> &values, const Array<double> &x) const
    {
        if (2*x[0] - x[1] > 0) values[0] = 1;
        //else if (x[0] < 0) values[0] = 2*x[0] + 1;
        else values[0] = 0;
    }
};

int main()
{
    //Create mesh
    std::array<Point, 2> a1 = {Point(-1, -1), Point(2, 1)};
    std::array<long unsigned int, 2> a2 = {100, 150};
    auto mesh = std::make_shared<Mesh>(RectangleMesh::create(a1, a2, CellType::Type::triangle));
    auto V = std::make_shared<Wave_equation::FunctionSpace>(mesh);

    //задание граничных условий
    auto u0 = std::make_shared<Boundary_Func>();
    auto boundary = std::make_shared<DirichletBoundary>();
    DirichletBC bc(V, u0, boundary);

    double dt = 0.0001;
    //создаём функции для расчёта
    auto u_pr = std::make_shared<Function>(V);
    auto u_prpr = std::make_shared<Function>(V);
    auto k = std::make_shared<Constant>(dt);
    auto ind = std::make_shared<Indicator>();

    //линейная и билинейная формы
    Wave_equation::BilinearForm a(V, V);
    Wave_equation::LinearForm L(V);
    auto u = std::make_shared<Function>(V);
    //auto g = std::make_shared<Constant>(0);
    a.dt = k;
    //L.dt = k;
    a.indicator = ind;
    //L.indicator = ind;
    for (int i = 0; i < 10000; i++)
    {
        u0->t = i * dt; //обозначение момента времени для граничной функции
        DirichletBC bc(V, u0, boundary); //создаём соответствующие граничные условия
        L.u_pr = u_pr;
        L.u_prpr = u_prpr;
        solve(a == L, *u, bc); //решаем

        std::cout << u0->t << '\n';
        if (i%5 == 0)
        {
        std::string fileName = "snapshots_plane_wave/Wave_equation-" + std::to_string(i/5) + ".pvd"; //сохраняем
        File file(fileName);
        file << *u;
        }
        *u_prpr = *u_pr; //переходим к новому моменту времени
        *u_pr = *u;
    }
    return 0;
}