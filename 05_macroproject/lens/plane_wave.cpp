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
        double c = 1.0/40;
        if (x[0] < c-1)
        {    if (t < 0.04) values[0] = 30 * sin(800 * t);}
        else values[0] = 0;   
    }
    double t;
};

class DirichletBoundary : public SubDomain //граница
{
    bool inside(const Array<double> &x, bool on_boundary) const
    {
        double c = 1.0 / 40;
        return ((x[0] < c - 1) or ((x[0] > -c-0.1) and (x[0] < c-0.1) and ((x[1] > 0.25) or (x[1] < -0.25))));
    }
};

class Indicator : public Expression //!!!Область с другим показателем преломления!!!
{
    void eval(Array<double> &values, const Array<double> &x) const
    {
        ///one
        ///if ((((2.7-x[0])*(2.7-x[0])) + x[1]*x[1] < 3*3) and ( x[0] <= 0)) values[0] = 1;
        ///if ((((2.7+x[0])*(2.7+x[0])) + x[1]*x[1] < 3*3) and ( x[0] >= 0)) values[0] = 1;
        ///two
        //if ((((9.9-x[0])*(9.9-x[0])) + x[1]*x[1] < 100) and ( x[0] <= 0)) values[0] = 1;
        //if ((((9.9+x[0])*(9.9+x[0])) + x[1]*x[1] < 100) and ( x[0] >= 0)) values[0] = 1;
        ///three
        if ((((1.4-x[0])*(1.4-x[0])) + x[1]*x[1] < 1.5*1.5) and ( x[0] <= 0)) values[0] = 1;
        if ((((1.4+x[0])*(1.4+x[0])) + x[1]*x[1] < 1.5*1.5) and ( x[0] >= 0)) values[0] = 1;
        ///four
        //if ((((1.4+x[0])*(1.4+x[0])) + x[1]*x[1] < 1.5*1.5) and ( x[0] >= 0)) values[0] = 1;
        ///five
        //if ((((2.7+x[0])*(2.7+x[0])) + x[1]*x[1] < 9) and ( x[0] >= 0)) values[0] = 1;
        else values[0] = 0;
    }
};

int main()
{
    //Create mesh
    std::array<Point, 2> a1 = {Point(-1, -1), Point(2, 1)};
    std::array<long unsigned int, 2> a2 = {200, 200};
    auto mesh = std::make_shared<Mesh>(RectangleMesh::create(a1, a2, CellType::Type::triangle));
    auto V = std::make_shared<Wave_equation::FunctionSpace>(mesh);

    //задание граничных условий
    auto u0 = std::make_shared<Boundary_Func>();
    auto boundary = std::make_shared<DirichletBoundary>();
    DirichletBC bc(V, u0, boundary);

    double dt = 0.00005;
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
    for (int i = 0; i < 5000; i++)
    {
        u0->t = i * dt; //обозначение момента времени для граничной функции
        DirichletBC bc(V, u0, boundary); //создаём соответствующие граничные условия
        L.u_pr = u_pr;
        L.u_prpr = u_prpr;
        solve(a == L, *u, bc); //решаем

        std::cout << u0->t << '\n';
        if (i%5 == 0)
        {
        std::string fileName = "snapshots_plane_wave3/Wave_equation-" + std::to_string(i/5) + ".pvd"; //сохраняем
        File file(fileName);
        file << *u;
        }
        *u_prpr = *u_pr; //переходим к новому моменту времени
        *u_pr = *u;
    }
    return 0;
}