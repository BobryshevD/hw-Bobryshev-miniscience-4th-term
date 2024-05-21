#include <dolfin.h>
#include <string>
#include "Initial_conditions.h"
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
        if (t < 0.04) values[0] = 30 * sin(1000 * t);
        else values[0] = 0;   
    }
    double t;
};

class DirichletBoundary : public SubDomain //граница
{
    bool wavering_zones (const Array<double> &x, double b, int N) const
    {
        bool indicator = false;
        for (int i = 0; i<N; i++)
            if ((x[1] > ((-(2*N+1)/2 + 2*i)*b)) && (x[1] < ((-(2*N+1)/2 + (2*i+1))*b))) indicator = true;
        return indicator;  
    }

    bool inside(const Array<double> &x, bool on_boundary) const
    {
        double c = 0.005;
        //return ((x[0] < c) && wavering_zones(x, 0.05, 5));
        return (x[0] < c);
    }
};

class Indicator : public Expression //экран
{
    void eval(Array<double> &values, const Array<double> &x) const
    {
        // double c = 0.02;
        // if (((x[1] > 0.1) or (x[1] < -0.1)) and ((x[0] > -0.5) and (x[0] < -0.2))) values[0] = c;
        // else if ((x[0] > -0.5) && ((x[1] > 0.85) || (x[1]) < -0.85)) values[0] = c;
        // else if (x[0] > 0.85) values[0] = c;
        // //else if (x[0] < 0) values[0] = 2*x[0] + 1;
        // else values[0] = 0;
        values[0] = 0;
    }
};

class Initial_condition: public Expression //начальные условия
{
    void eval(Array<double> &values, const Array<double> &x) const
    {
        values[0] = 0;
    }
};

int main()
{
    //Create mesh
    // std::array<Point, 2> a1 = {Point(-1, -1), Point(1, 1)};
    // std::array<long unsigned int, 2> a2 = {100, 100};
    auto mesh = std::make_shared<Mesh>("test_mesh.xml");
    
    
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
    auto ic = std::make_shared<Initial_condition>();

    //линейная и билинейная формы
    Wave_equation::BilinearForm a(V, V);
    Wave_equation::LinearForm L(V);

    Initial_conditions::BilinearForm a_ic(V, V);
    Initial_conditions::LinearForm L_ic(V);
    L_ic.ic = ic;

    solve(a_ic == L_ic, *u_pr);
    solve(a_ic == L_ic, *u_prpr);

    auto u = std::make_shared<Function>(V);
    //auto g = std::make_shared<Constant>(0);
    a.dt = k;
    L.dt = k;
    a.indicator = ind;
    L.indicator = ind;
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